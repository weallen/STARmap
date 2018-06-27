from scipy.spatial import ConvexHull
from skimage.transform import downscale_local_mean
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from __future__ import print_function
from skimage.measure import regionprops

def GetQHulls(labels):
    Nlabels = labels.max()
    hulls = []
    coords = []
    num_cells = 0
    for i, region in enumerate(regionprops(labels)):    
        print(i,"/",Nlabels)
        curr_coords = region.coords #np.argwhere(labels == i)
        # size threshold of > 100 pixels and < 100000
        if curr_coords.shape[0] < 100000 and region.area > 100:
            num_cells += 1
            hulls.append(ConvexHull(curr_coords))
            coords.append(curr_coords)
    print("Used %d / %d" % (num_cells, Nlabels))
    return hulls, coords


def hull_to_polygon(hull):    
    cent = np.mean(hull.points, 0)
    pts = []
    for pt in hull.points[hull.simplices]:
        pts.append(pt[0].tolist())
        pts.append(pt[1].tolist())
    pts.sort(key=lambda p: np.arctan2(p[1] - cent[1],
                                    p[0] - cent[0]))
    pts = pts[0::2]  # Deleting duplicates
    pts.insert(len(pts), pts[0])
    k =1.1
    poly = Polygon(k*(np.array(pts)- cent) + cent,edgecolor='k', linewidth=1)
    #poly.set_capstyle('round')
    return poly

def plot_poly_cells_expression(nissl, hulls, expr, cmap, good_cells=None,width=2, height=9,figscale=10,alpha=1):
    figscale = 10
    plt.figure(figsize=(figscale*width/float(height),figscale))
    polys = [hull_to_polygon(h) for h in hulls]
    if good_cells is not None:
        polys = [p for i,p in enumerate(polys) if i in good_cells]
    p = PatchCollection(polys,alpha=alpha, cmap=cmap,linewidths=0)
    p.set_array(expr)
    p.set_clim(vmin=0, vmax=max(colors))        
    plt.gca().add_collection(p)
    plt.axis('off')

def plot_poly_cells_cluster(nissl, hulls, colors, cmap, good_cells=None,width=2, height=9,figscale=10, rescale_colors=False,alpha=1,vmin=None,vmax=None):
    figscale = 10
    plt.figure(figsize=(figscale*width/float(height),figscale))
    polys = [hull_to_polygon(h) for h in hulls]
    if good_cells is not None:
        polys = [p for i,p in enumerate(polys) if i in good_cells]
    p = PatchCollection(polys,alpha=alpha, cmap=cmap,edgecolor='k', linewidth=0.5)
    if vmin or vmax is not None:
        p.set_array(colors)
        p.set_clim(vmin=vmin,vmax=vmax)
    else:
        if rescale_colors:
            p.set_array(colors+1)
            p.set_clim(vmin=0, vmax=max(colors+1))
        else:
            p.set_array(colors)
            p.set_clim(vmin=0, vmax=max(colors))        
    nissl = (nissl > 0).astype(np.int)
    plt.imshow(nissl.T,cmap=plt.cm.gray_r,alpha=0.15)
    plt.gca().add_collection(p)
    plt.axis('off')
    return polys

def plot_cells_cluster(nissl, coords, good_cells, colors, cmap, width=2, height=9,figscale=100, vmin=None,vmax=None):
    figscale = 10
    plt.figure(figsize=(figscale*width/float(height),figscale))
    img = -1*np.ones_like(nissl)
    curr_coords = [coords[k] for k in range(len(coords)) if k in good_cells]
    for i,c in enumerate(curr_coords):
        for k in c:
            if k[0] < img.shape[0] and k[1] < img.shape[1]:
                img[k[0],k[1]] = colors[i]
    plt.imshow(img.T,cmap=cmap,vmin=-1,vmax=colors.max())
    plt.axis('off')
                      
def get_cells_and_clusts_for_experiment(analysis_obj, expt_id):
    good_cells = analysis_obj._meta.index[(analysis_obj._meta["orig_ident"]==expt_id)].values
    colors = analysis_obj._clusts[analysis_obj._meta["orig_ident"]==expt_id]
    return good_cells, colors
