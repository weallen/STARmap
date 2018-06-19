#
# CODE FOR CONVERTING PROCESSED READS + NISSL STAIN PER CELL QUANTIFICATION
#
#
import matplotlib
import tifffile
import os
import re
import sys
import numpy as np
import matplotlib.pyplot as plt
from skimage.feature import peak_local_max
from scipy.spatial.distance import cdist
from scipy.ndimage.filters import gaussian_laplace
from skimage.transform import downscale_local_mean, SimilarityTransform, warp
from skimage.filters import laplace, gaussian
from skimage.morphology import binary_erosion
import SimpleITK as sitk
from scipy.spatial import cKDTree
from pandas import DataFrame
from joblib import *
import collections
import sys
from scipy.io import loadmat
from coding import *
from scipy.spatial import ConvexHull
import seaborn as sns
import xml.etree.ElementTree as ET

def ParseCellCounter(path):
    """
    Parse cell locations exported from Fiji CellCounter.
    Used to export manually selected cell locations from clicking on DAPI+ nuclei.
    """

    tree = ET.parse(path)
    root = tree.getroot()[1]
    vals = []
    for i, child in enumerate(root[1].findall("Marker")):
        x = int(child[0].text)
        y = int(child[1].text)
        vals.append([x,y])
    return np.array(vals)

def LoadIlastikImage(fpath):
    """Loads an Ilastik exported image, filters with gaussian, then thresholds.
    
    Arguments:
        fpath: path of image
    
    Returns:
        img: binary image where predicted > threshold
    """

    from skimage.filters import gaussian
    img = (tifffile.imread(fpath)-1.)*255.
    img = gaussian(img, sigma=2)
    return img>250

def LoadNisslData(dirname,fname="nissl_maxproj_resized.tif"):
    """Load Nissl data from directory containing nissl subdirectory.
    Assumes directory structure is:
        fpath/nissl/nissl_maxproj_resized.tif 
    Arguments:
        fpath: path to 
    
    Returns:
        nissl: nissl image from directory
    """

    nissl = tifffile.imread(os.path.join(dirname, "nissl", fname))
    return nissl

def LoadCellPoints(fpath):
    S = loadmat(os.path.join(fpath, "output", "cellLocs.mat"))
    return np.round(S["cellLocs"])

def LoadReadPos(fpath):
    S = loadmat(os.path.join(fpath, "output", "goodPoints.mat"))
    bases = [str(i[0]) for i in S["goodBases"][0]]
    points = S["goodPoints"][:,:2]
    temp = np.zeros(points.shape)
    temp[:,0] = np.round(points[:,1]-1)
    temp[:,1] = np.round(points[:,0]-1)
    return bases, temp

def LoadGenes(fpath):
    genes2seq = {}
    seq2genes = {}
    with open(os.path.join(fpath, "genes.csv")) as f:
        for l in f:
            fields = l.rstrip().split(",")
            genes2seq[fields[0]] = "".join([str(s+1) for s in EncodeSOLID(fields[1][::-1])])
            seq2genes[genes2seq[fields[0]]] = fields[0]
    return genes2seq, seq2genes

def SegmentNisslData(fpath, nissl, cell_locs):
    # uses watershed to try this
    from skimage.morphology import watershed, binary_dilation
    from scipy import ndimage as ndi
    from skimage.morphology import watershed
    import cv2
    from skimage.morphology import disk

    from skimage.feature import peak_local_max
    blurred_nissl_seg = gaussian(nissl.astype(np.float),10) > 50
    print("Dilating")
    blurred_nissl_seg = binary_dilation(blurred_nissl_seg, selem=disk(10))
    print("Distance transform")
    markers = np.zeros(blurred_nissl_seg.shape, dtype=np.uint8)
    for i in range(cell_locs.shape[0]):
        y,x = cell_locs[i,:]
        if x < blurred_nissl_seg.shape[0] and y < blurred_nissl_seg.shape[1]:
            markers[x-1,y-1] = 1
    markers = ndi.label(markers)[0]
    print("Watershed")
    labels = watershed(blurred_nissl_seg, markers, mask=blurred_nissl_seg)
    labels_line = watershed(blurred_nissl_seg, markers, mask=blurred_nissl_seg,watershed_line=True)
    print("Labeled %d cells" % labels.max())
    tifffile.imsave(os.path.join(fpath, "output", "labeled_cells_line.tif"),labels_line.astype(np.uint16))
    tifffile.imsave(os.path.join(fpath, "output", "labeled_cells.tif"), labels.astype(np.uint16))
    return labels

def AssignReadsToCells(fpath, labels, good_spots, genes):
    Nlabels = labels.max()
    # make matrix of  XY coordinates of label pixels and corresponding labels
    Npixels = len(np.where(labels > 0)[0])
    coords = []
    cell_ids = []
    print("Grabbing coordinates of cells")
    num_cells = 0
    for i in range(Nlabels): # skip label 0 (background)
        curr_coords = np.argwhere(labels == i)
        if curr_coords.shape[0] < 100000 and curr_coords.shape[0] > 1000:
            coords.append(curr_coords)
            cell_ids.append(np.repeat(i, curr_coords.shape[0]))
            num_cells += 1
        else:
            coords.append(np.array([[],[]]).T)
    print("Using %d out of %d cells" % (num_cells, Nlabels))
    coords_list = coords
    coords = np.vstack(coords)
    cell_ids = np.concatenate(cell_ids)
    print("Building KD tree of cell coords")
    label_kd = cKDTree(coords)
    print("Assigning reads to cells")
    #print query_results[:10]
    cell_assignments = np.array([cell_ids[label_kd.query(p)[1]] for p in good_spots])  #
    return cell_assignments, coords_list


def GetQHulls(fpath, labels):
    Nlabels = labels.max()
    hulls = []
    coords = []
    num_cells = 0
    for i in range(Nlabels): # skip label 0 (background)
        curr_coords = np.argwhere(labels == i)
        if curr_coords.shape[0] < 100000 and curr_coords.shape[0] > 1000:
            num_cells += 1
            hulls.append(ConvexHull(curr_coords))
            coords.append(curr_coords)
    return hulls, coords

def AssignReadsToCellsQHull(fpath, labels, good_spots):
    from matplotlib.path import Path
    Nlabels = labels.max()
    # make matrix of  XY coordinates of label pixels and corresponding labels
    print("Grabbing coordinates of cells")
    num_cells = 0
    hulls = []
    coords = []
    for i in range(Nlabels): # skip label 0 (background)
        curr_coords = np.argwhere(labels == i)
        if curr_coords.shape[0] < 100000 and curr_coords.shape[0] > 1000:
            num_cells += 1
            hulls.append(ConvexHull(curr_coords))
            coords.append(curr_coords)
#    for i in range(len(hulls)):
    print("Assigning reads to cells")
    point_assignments = []

    for i, h in enumerate(hulls):
        p = Path(h.points[h.vertices])
        point_assignments.append(np.argwhere(p.contains_points(good_spots)).flatten()) 
    return hulls, point_assignments, coords


def ConvertReadAssignmentsQHull(fpath, point_assignments, bases, seqs2genes):
    outdir = os.path.join(fpath, "output", "singlecell")
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    gene_seqs = seqs2genes.keys()
    Ncells = len(point_assignments)
    cell_by_barcode = np.zeros((Ncells,len(gene_seqs)))
    gene_seq_to_index = {} # map from sequence to index into matrix

    for i,k in enumerate(gene_seqs):
        gene_seq_to_index[k] = i

    print(gene_seq_to_index.keys())
    print("Counting reads")
    total_read_count = 0
    for i in range(Ncells):
        if i % 50 == 0:
            print("Cell %d" % i)
        assigned_barcodes = point_assignments[i] # which peaks are assigned to that cell
        for j in assigned_barcodes: # which actual colorseq those correspond t
            b = bases[j]
            if b in gene_seq_to_index:
                cell_by_barcode[i,gene_seq_to_index[b]] += 1
                total_read_count += 1

    #print "%f percent [%d out of %d] reads were assigned to cells" % (total_read_count/Ngood, total_read_count, Ngood)
    np.save(os.path.join(outdir, "cell_barcode_count.npy"), cell_by_barcode)
    np.savetxt(os.path.join(outdir, "cell_barcode_count.csv"), cell_by_barcode.astype(np.int), delimiter=',', fmt="%d")
    f = open(os.path.join(outdir, "cell_barcode_names.csv"),'w')
    for i,k in enumerate(gene_seqs):
        f.write("%d,%s,%s\n" % (i, k, seqs2genes[k]))
    f.close()
    return cell_by_barcode

def ConvertReadAssignments(fpath, good_spots, Nlabels, cell_assignments, bases, seqs2genes):
    outdir = os.path.join(fpath, "output", "singlecell")
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    gene_seqs = seqs2genes.keys()
    #Nlabels = cell_assignments.flatten().max()
    cell_by_barcode = np.zeros((Nlabels,len(gene_seqs)))
    gene_seq_to_index = {}

    for i,k in enumerate(gene_seqs):
        gene_seq_to_index[k] = i

    print(gene_seq_to_index.keys())
    print("Counting reads")
    total_read_count = 0
    for i in range(Nlabels):
        if i % 50 == 0:
            print("Cell %d" % i)
        assigned_barcodes = np.where(cell_assignments==i)[0] # which peaks are assigned to that cell
        for j in assigned_barcodes: # which actual colorseq those correspond t
            b = bases[j]
            cell_by_barcode[i,gene_seq_to_index[b]] += 1
            total_read_count += 1

    Ngood = float(good_spots.shape[0])
    print("%f percent [%d out of %d] reads were assigned to cells" % (total_read_count/Ngood, total_read_count, Ngood))
    np.save(os.path.join(outdir, "cell_barcode_count.npy"), cell_by_barcode)
    np.savetxt(os.path.join(outdir, "cell_barcode_count.csv"), cell_by_barcode.astype(np.int), delimiter=',', fmt="%d")
    f = open(os.path.join(outdir, "cell_barcode_names.csv"),'w')
    for i,k in enumerate(gene_seqs):
        f.write("%d,%s,%s\n" % (i, k, seqs2genes[k]))
    f.close()
    return cell_by_barcode

def SaveExpressionImages(fpath, d, labels, hulls):
    from scipy.misc import imresize

    outdir = os.path.join(fpath, "output", "singlecell")
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    gene_names = d.columns
    for g in gene_names:
        img = MakeExpressionImage(g, d, labels, hulls)
        plt.figure(figsize=(20,20))
        plt.imshow(imresize(img, 0.25), cmap=plt.cm.jet)
        plt.axis('off')
        plt.savefig(os.path.join(outdir, g+"_cells.png"))
        plt.close()

def MakeExpressionImage(gene_name, d, labels, hulls):    
    Nlabels = len(hulls)
    expr = d[gene_name]
    expr_img = np.zeros_like(labels)
    for i in range(Nlabels):
        p = hulls[i].points.astype(np.int)
        expr_img[p[:,0], p[:,1]] = expr[i]
    return expr_img


def PlotCellNumbers(fpath, labels):
    from skimage.measure import regionprops
    outdir = os.path.join(fpath, "output")
    plt.figure(figsize=(20,10))
    plt.imshow(labels,cmap=plt.cm.jet)
    #for i in range(cell_locs.shape[0]):
    for i, region in enumerate(regionprops(labels)):
        plt.text(region.centroid[1], region.centroid[0], str(i), fontsize=7, color='w')
    plt.savefig(os.path.join(outdir, "cell_nums.png"))

def PlotClusters(fpath, labels, hulls, ident_name, outname, cmap=None):
    import matplotlib.patches as mpatches

    outpath = os.path.join(fpath, "output")
    # load cluster labels
    num2ident = {}
    max_ident = 0
    with open(os.path.join(outpath, ident_name)) as f:
        for i,l in enumerate(f):
            if i > 0:
                name, ident = l.rstrip().split(",")
                cell_num = int(name.split("_")[1])+1
                num2ident[cell_num] = int(ident)
                if int(ident) > max_ident:
                    max_ident = int(ident)
    cluster_img = np.zeros_like(labels)
    for k,v in num2ident.items():
        p = hulls[k-1].points.astype(np.int)
        cluster_img[p[:,0], p[:,1]] = v+1
    plt.figure(figsize=(20,10))
    if cmap is None:
        cmap = plt.cm.OrRd
    values = range(cluster_img.max()+1)#[0,1,2,3,4,5]
    im = plt.imshow(cluster_img, cmap=cmap,vmin=0,vmax=max_ident+1)
    colors = [ im.cmap(im.norm(value)) for value in values]
    patches = [ mpatches.Patch(color=colors[i], label="Cluster {l}".format(l=values[i]-1) ) for i in range(1,len(values)) ]
    plt.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0. )
    plt.axis('off')
    plt.savefig(os.path.join(outpath, outname), transparent=True)

def main():
    # input is list of absolute paths to folders for processing
    paths = sys.argv[1:]

    dropbox_path = "/home/dlab/Dropbox/snailseq/data"
    use_ilastik = False
    assign_reads = True
    for i in range(len(paths)):
        fpath = paths[i]
        print("\t Processing %s" % fpath)
        path = paths[i]
        genes2seq, seq2genes = LoadGenes(fpath)
        bases, points = LoadReadPos(fpath)

        cell_locs = ParseCellCounter(os.path.join(dropbox_path, "dapi", "CellCounter_" + path + ".xml"))
        nissl = LoadNisslData(fpath)
        labels = SegmentNisslData(fpath, nissl, cell_locs.astype(np.int))
        PlotCellNumbers(fpath, labels)
        plt.imsave(fname=os.path.join(fpath, "output", "labels.tif"), arr=labels)

        # point points + segmentation
        plt.figure(figsize=(80,40))
        plt.plot(points[:,1], points[:,0],'r.',markersize=0.5)
        plt.imshow(labels,cmap=plt.cm.gray)
        plt.axis('off')
        points_seg_path = os.path.join(fpath, "output", "points_seg.png")
        print("Saving %s" % points_seg_path)
        plt.savefig(points_seg_path)
        if assign_reads:
            hulls, point_assignments, coords = AssignReadsToCellsQHull(fpath, labels, points)
            cell_by_barcode = ConvertReadAssignmentsQHull(fpath, point_assignments, bases, seq2genes)
            print("Saving out")

            pts = [s.points for s in hulls]
            clusters_path = os.path.join(dropbox_path, "clusters", paths[i])
            if not os.path.exists(clusters_path):
                os.mkdir(clusters_path)
            np.savez(os.path.join(clusters_path, "labels.npz"), labels=labels)


if __name__ == "__main__":
    main()
