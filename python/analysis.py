
from __future__ import print_function
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA,FactorAnalysis,NMF
from sklearn.manifold import TSNE
from sklearn.cluster import SpectralClustering, DBSCAN
from sklearn.preprocessing import scale, MinMaxScaler
from sklearn.neighbors import NearestNeighbors
from scipy.stats import ttest_ind, norm, ranksums
from scipy.stats.mstats import gmean

from skimage.transform import downscale_local_mean

import os
from sklearn.linear_model import LinearRegression 
import umap

import statsmodels.formula.api as smf
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

def load_data(data_dir, prefix="Cell"):
        #expr = pd.read_csv(os.path.join(data_dir, "data_table.csv"), index_col=0)
    expr = pd.read_csv(os.path.join(data_dir, "cell_barcode_count.csv"), header=None)
    gene_names = pd.read_csv(os.path.join(data_dir, "cell_barcode_names.csv"),header=None)
    rownames = ["%s_%05d"% (prefix,i) for i in range(expr.shape[0])]
    names = gene_names[2]
    names.name = "Gene"
    return pd.DataFrame(data=expr.values, columns=names, index=rownames)

class STARmapAnalysis(object):    
        
    def __init__(self):
        self._raw_data = None # raw data
        self._data = None # data that has been normalized (log + total count)
        self._scaled = None # scaled data
        self._pca = None # pca for cells
        self._tsne = None # tsne for cells
        self._clusts = None # per cell clustering
        self._meta = None # per cell metadata
        self._nexpt = 0

        self._nissl = None
        self._hulls = None

    # basic usage:
    # S = STARmapAnalysis()
    # S.load_data(data_dir)
    # S.scale()
    # S.pca()
    # S.tsne()
    # S.cluster()
    # DATA LOADING FUNCTIONS
    #
    
        
#        self.add_data()
        # rows are cells

    def get_cells_by_experiment(self, idx, use_genes=None,use_scaled=False):
        expt_idx = np.argwhere(self._meta["orig_ident"]==idx).flatten()
        if use_scaled:
            data = self._raw_data
        else:
            data = self._data
        if use_genes:
            return data.iloc[expt_idx,:].loc[:,use_genes]
        else:
            return data.iloc[expt_idx,:]

    def subset_cells(self, cell_names):
        pass

    def add_data(self, data, group=None,tpm_norm=False,use_genes=None):
        """
        Add data to data matrix, keeping track of source.
        :input min_genes: include cells expressing >= min_genes
        :input min_cells: include genes expressed in >= min_cells 
        """
        
        
        reads_per_cell = data.sum(axis=1)
        genes_per_cell = (data > 0).sum(axis=1)
        cells_per_gene = (data > 0).sum(axis=0)

        if tpm_norm:
            data = data / (data.sum().sum()/1e6)
        if use_genes:
            data = data[use_genes]
        if group is not None: 
                meta = pd.DataFrame(np.vstack((np.ones((len(reads_per_cell),),dtype=np.int)*self._nexpt, reads_per_cell, genes_per_cell, np.repeat(group, len(reads_per_cell)))).T, \
                                      columns=["orig_ident", "reads_per_cell", "genes_per_cell", "group"])
        else:
            meta = pd.DataFrame(np.vstack((np.ones((len(reads_per_cell),),dtype=np.int)*self._nexpt, reads_per_cell, genes_per_cell)).T, \
                                      columns=["orig_ident", "reads_per_cell", "genes_per_cell"])
        
        if self._nexpt == 0:
            self._raw_data = data
            self._meta = meta
        else: # add data to existing dataframe
            self._raw_data = self._data.append(data)
            self._meta = self._meta.append(meta)
        self._data = self._raw_data
        self._nexpt += 1
        self._ncell, self._ngene = self._data.shape
            

    
    #
    # PREPROCESSING FUNCTIONS
    #
    def filter_cells_by_expression(self, min_genes = 10, min_cells = 10):
        good_cells = (self._raw_data.values > 0).sum(axis=1)>min_genes 
        good_genes = (self._raw_data.values > 0).sum(axis=0)>min_cells
        self._raw_data = self._raw_data.loc[good_cells,:]
        self._raw_data = self._raw_data.loc[:,good_genes]
        self._data = self._data.loc[good_cells,:]
        self._data = self._data.loc[:,good_genes]
        self._meta = self._meta.loc[good_cells]

    def filter_cells_by_feature(self, feature_name, low_thresh, high_thresh):
        m = self._meta[feature_name].values
        to_keep = np.logical_and(m > low_thresh, m <= high_thresh)
        self._raw_data = self._raw_data.loc[to_keep,:]
        self._data = self._data.loc[to_keep,:]
        if self._scaled is not None:
            self._scaled = self._scaled.loc[to_keep,:]
        self._meta = self._meta.loc[to_keep,:]
        self._ncell, self._ngene = self._data.shape

    def median_normalize(self):
        data = self._data.values.copy().astype(np.float)
        for i in range(data.shape[1]):
            nonzero = np.argwhere(data[:,i] != 0).flatten()
            data[nonzero,i] =  data[nonzero,i]/gmean(data[nonzero,i])
        norm_data = np.zeros_like(data)
        for i in range(data.shape[0]):
            temp = data[i,:]
            norm_data[i,:] = np.median(temp[temp!=0])
        for i in range(data.shape[0]):
            data[i,:] /= norm_data[i,:]
        self._data = pd.DataFrame(data=data, columns=self._data.columns,index=self._data.index)
        #self._data.values = data

    def normalize(self, norm_type="none" , use_genes=None, scale_factor=10000):
        """ Normalize  for clustering and TPM normalize raw_data """
        if use_genes:
            data = self._data.loc[:,use_genes]
        else:
            data = self._data
        median_transcripts = np.median(self._raw_data.sum(axis=1))
        for i in range(self._ncell):
            if norm_type == "abs":
                self._data.iloc[i,:] = np.log1p((self._data.iloc[i,:] / data[i,:].sum()) * scale_factor)
            elif norm_type == "median":
                # normalize to the median transcripts per cell, across all cells
                self._data.iloc[i,:] = np.log1p(median_transcripts * (self._data.iloc[i,:] / data.iloc[i,:].sum()))
            elif "none":
                self._data.iloc[i,:] = np.log1p(self._data.iloc[i,:])

    def scale(self, vars_to_regress=None, model_type="none", do_trim=False, do_scale=True, do_center=True, scale_max = 10):
        """
        Regress out reads per cell and identity
        """
        scaled = np.zeros((self._ncell, self._ngene))
        reads_per_cell = self._meta["reads_per_cell"]
        genes_per_cell = self._meta["genes_per_cell"]
        ident = self._meta["orig_ident"]
        group = self._meta["group"]
        if model_type is "none":
            scaled = self._data.values.copy()
        else:
            for i in range(self._ngene):
                expr = self._data.iloc[:,i]
                d = pd.DataFrame(np.array((expr.astype(np.float),reads_per_cell, genes_per_cell, ident, group)).T,columns=["expr", "reads_per_cell", "genes_per_cell", "orig_ident", "group"])            
                if model_type is "linear":
                    results = smf.ols('expr ~ orig_ident + reads_per_cell + group', data=d).fit()
                    scaled[:,i] = results.resid
                elif model_type is "poisson":
                    results = smf.glm('expr ~ reads_per_cell + orig_ident + group', data=d,family=sm.families.Poisson()).fit()
                    #results = smf.glm('expr ~ orig_ident', data=d,family=sm.families.NegativeBinomial()).fit()
                    scaled[:,i] = results.resid_pearson

        self._scaled = pd.DataFrame(scaled, columns=self._data.columns, index=self._data.index)
        if do_trim:
            x = self._scaled.mean()
            y = self._scaled.var()/x
            plt.plot(x,y,'.')
            good_genes = np.array(np.logical_and(y.values>1, x.values>0.1))
            self._scaled = self._scaled.iloc[:,good_genes]

        if do_center or do_scale:
            for i in range(self._scaled.shape[1]):
                temp = self._scaled.iloc[:,i].values
                temp = scale(temp, with_mean=do_center, with_std=do_scale)
                temp[temp > scale_max] = scale_max
                self._scaled.iloc[:,i] = temp 
    #
    # DATA ACCESS FUNCTIONS
    #
    def get_gene_names(self):
        return self._raw_data.columns

    def get_metadata_names(self):
        return self._meta.columns

    def get_metadata(self, m):
        return self._meta[m]
    
    def get_expr_for_gene(self, gene_name, scaled=True):
        if scaled:
            return np.log1p(self._raw_data[gene_name])
        else:
            return self._raw_data[gene_name]
    
    def get_expr_for_cell(self, i, scaled=True):
        if scaled:
            return np.log1p(self._raw_data.iloc[i,:])
        else:
            return self._raw_data.iloc[i,:]    

    def get_mean_expression_across_clusters(self, scaled=True):
        expr = [self.get_mean_expr_for_cluster(i) for i in range(max(self._clusts))]
        return pd.concat(expr,axis=1).transpose()

    def get_mean_expr_for_cluster(self, clust_id, scaled=True):
        data = self.get_cells_by_cluster(clust_id, use_raw=True)
        if scaled:
            return np.log1p(data)
        else:
            return data
    
    def get_cells_by_cluster(self, clust_id, use_raw=False):
        cells = self._clusts == clust_id
        if use_raw:
            return self._raw_data.iloc[cells,:]
        else:
            return self._data.iloc[cells,:]

    def get_metadata_by_cluster(self, clust_id):
        cells = self._clusts == clust_id
        return self._meta.iloc[cells,:]

    def get_clusts_for_experiment(self, expt_id):
        return self._clusts[self._meta["orig_ident"]==expt_id]
    #
    # BASIC QC FUNCTIONS
    # 

    def merge_multiple_clusters(self, clusts):
        # merge each cluster into first in list
        for idx,c in enumerate(clusts[1:]):
            self._clusts[self._clusts==c] = clusts[0]
        temp = self._clusts.copy()
        # relabel clusters to be contiguous
        for idx, c in enumerate(np.unique(self._clusts)):
            temp[self._clusts==c] = idx
        self._clusts = temp

    def merge_clusters(self, clust0, clust1):
        self._clusts[self._clusts==clust1] = clust0
        temp = self._clusts.copy()
        for idx,c in enumerate(np.unique(self._clusts)):
            temp[self._clusts==c] = idx
        self._clusts = temp

    def subset_by_cellidx(self, cell_ids):
        S = SnailSeqAnalysis()
        S._raw_data = self._raw_data.iloc[cell_ids,:] # raw data
        S._data = self._data.iloc[cell_ids,:] # data that has been normalized (log + total count)
        S._scaled = self._scaled.iloc[cell_ids,:] # scaled data
        S._meta = self._meta.iloc[cell_ids] # per cell metadata
        S._nexpt = self._nexpt
        S._ncell, S._ngene = S._data.shape
        return S

    def subset_by_cluster(self, clust):
        cell_ids = np.argwhere(np.array(self._clusts) == clust).flatten()
        S = self.subset_by_cellidx(cell_ids)
        return S

    def plot_stats_per_cell(self,gene_thresh=0):
        reads_per_cell = self._meta["reads_per_cell"]
        genes_per_cell = self._meta["genes_per_cell"]

        plt.subplot(2,3,1)
        plt.hist(reads_per_cell,10, color='k')
        plt.ylabel('# cells')
        plt.xlabel('# reads')

        plt.subplot(2,3,2)
        plt.hist(genes_per_cell,10,color='k')
        plt.ylabel('# cells')
        plt.xlabel('# genes') 

        plt.subplot(2,3,3)
        plt.title('R=%f' % np.corrcoef(reads_per_cell.T, genes_per_cell)[0,1])
        plt.scatter(reads_per_cell, genes_per_cell,marker='.',s=30, c=self._meta["orig_ident"],cmap=plt.cm.jet,lw=0)
        plt.xlabel("Reads per cell")
        plt.ylabel("Genes per cell")




    # DIMENSIONALITY REDUCTION
    #
    def pca(self, n_components=10, use_genes=None, use_corr=False):
        if use_genes is not None:
            d = self._scaled.loc[:,use_genes].dropna(axis=1)
        else:
            d = self._scaled
        if use_corr:
            d = np.corrcoef(d)
        if use_genes is None:
            self._pca = PCA(n_components=n_components).fit(d)
            self._transformed_pca = self._pca.transform(d)
        else:
            self._pca = PCA(n_components=n_components).fit(d)
            self._transformed_pca = self._pca.transform(d)

    def umap(self, max_pc=5, n_neighbors=10, min_dist=0.3, metric="euclidean"):
        self._tsne = umap.UMAP(n_neighbors=n_neighbors, min_dist=min_dist, metric="euclidean").fit_transform(self._transformed_pca[:,:max_pc])

    def tsne(self, max_pc=5, perplexity=20):
        self._tsne = TSNE(n_components=3, random_state=1, perplexity=perplexity).fit_transform(self._transformed_pca[:,:max_pc])

    def plot_tsne(self,cmap=None, s=10, renumber_clusts=False):
       
        if self._clusts is None:
            plt.plot(self._tsne[:,0], self._tsne[:,1],'o')
        else:
            if cmap is None:
                cmap = plt.cm.get_cmap('jet', len(np.unique(self._clusts)))
            if not renumber_clusts:
                plt.scatter(self._tsne[:,0], self._tsne[:,1],c=self._clusts,s=s,cmap=cmap,lw=0)
            else:
                plt.scatter(self._tsne[:,0], self._tsne[:,1],c=self._clusts+1,s=s,cmap=cmap,vmin=0, vmax=self._clusts.max()+1,lw=0)
        
    def plot_pca(self, dim0=0, dim1=1,s=10,cmap=plt.cm.jet):
        if self._clusts is None:
            plt.scatter(self._transformed_pca[:,dim0], self._transformed_pca[:,dim1],c=self.get_metadata("orig_ident"), cmap=cmap,s=s,lw=0)
        else:
            plt.scatter(self._transformed_pca[:,dim0], self._transformed_pca[:,dim1],c=self._clusts,cmap=cmap,s=s,lw=0)

    #
    # CLUSTERING
    #
    def cluster_hdbscan(self,max_pc=5):
        import hdbscan
        #clusterer = hdbscan.RobustSingleLinkage(cut=0.125, k=7)
        clusterer = hdbscan.HDBSCAN(min_cluster_size=10,alpha=2.)
        self._clusts = np.array(clusterer.fit_predict(self._tsne))

    def cluster_dbscan(self, eps=0.5):
        self._clusts = np.array(DBSCAN(eps=eps).fit_predict(self._tsne))
    
    def cluster_snn(self, max_pc=5, k=30):
        import igraph as ig
        data = self._transformed_pca[:,:max_pc]
        nbrs = NearestNeighbors(n_neighbors=k, algorithm='kd_tree').fit(data)
        neighbor_graph = nbrs.kneighbors_graph(data)
        g = ig.Graph()
        g = ig.GraphBase.Adjacency(neighbor_graph.toarray().tolist(), mode=ig.ADJ_UNDIRECTED)
        sim = np.array(g.similarity_jaccard())
        g = ig.GraphBase.Weighted_Adjacency(sim.tolist(), mode=ig.ADJ_UNDIRECTED)
        self._clusts = np.array(g.community_multilevel(weights="weight", return_levels=False))


    def cluster_gmm(self, n_clusts=5,max_pc=5):

        from sklearn.mixture import GaussianMixture
        model = GaussianMixture(n_components=n_clusts).fit(self._transformed_pca[:,:max_pc])
        self._clusts = np.array(model.predict(self._transformed_pca[:,:max_pc]))

    def load_clustering_data(self, path,encoding=None):
        """
        Load clusters and dimensionality reduction from .npy file
        """
        if encoding is None:
            data = np.load(path)
        else:
            data = np.load(path,encoding=encoding)
        self._tsne = data["tsne"]
        self._pca = data["pca"]
        self._clusts = data["clusts"]

    def save_clusters(self, experiment_id, outpath):
        f = open(outpath,'w')
        f.write("names,idents\n")
        cell_names = list(self.get_cells_by_experiment(experiment_id).index)
        clust_idx = np.argwhere(self._meta["orig_ident"]==experiment_id).flatten()
        for i,c in enumerate(self._clusts[clust_idx]):
            f.write("%s,%d\n" % (cell_names[i],c))

    def plot_explained_variance(self):
        if self._pca is not None:
            plt.plot(self._pca.explained_variance_ratio_, 'ko-')
            
    def find_all_markers(self, test="bimod", use_genes=None, only_pos=True, log_fc_thresh=0.25, min_pct=0.1, fdr_correct=True):
        dfs = []
        for i in np.unique(self._clusts):
            if i >= 0:
                markers = self.find_markers(i, test=test, use_genes=use_genes, only_pos=only_pos, log_fc_thresh=log_fc_thresh, min_pct=min_pct)
                markers['cluster'] = i 
                dfs.append(markers.sort_values(["pval", "log_fc"]))
        dfs = pd.concat(dfs)
        if fdr_correct:
            _, qvals, _, _ = multipletests(dfs["pval"], method="fdr_bh")
            dfs["pval"] = qvals
        return dfs


    def find_markers(self, clust0, clust1=None, test="bimod", use_genes=None, only_pos=True, log_fc_thresh=0.25, min_pct=0.1):
        if use_genes is None:
            curr_data = self._data
        else:
            curr_data  = self._data.loc[:,use_genes]
        cells1 = np.argwhere(self._clusts == clust0).flatten()
        if clust1 is not None:
            cells2 = np.argwhere(self._clusts == clust1).flatten()
        else:
            cells2 = np.argwhere(self._clusts != clust0).flatten()

        # select genes based on being expressed in a minimum fraction  of cells
        fraction1 = (curr_data.iloc[cells1,:]>0).sum(axis=0)/float(len(cells1))
        fraction2 = (curr_data.iloc[cells2,:]>0).sum(axis=0)/float(len(cells2))
        good_frac = np.logical_or(fraction1>min_pct, fraction2>min_pct)
        fraction1 = fraction1[good_frac]
        fraction2 = fraction2[good_frac]
        curr_data = curr_data.loc[:,good_frac]

        # select genes based on FC
        log_fc = np.array([np.log1p(np.expm1(curr_data.iloc[cells1,i].values).mean()) - np.log1p(np.expm1(curr_data.iloc[cells2,i].values).mean()) for i in range(curr_data.shape[1])])
        if only_pos:
            good_fc = log_fc>log_fc_thresh
        else:
            good_fc = np.abs(log_fc)>log_fc_thresh
        curr_data = curr_data.iloc[:,good_fc]
        log_fc = log_fc[good_fc]
        fraction1 = fraction1[good_fc]
        fraction2 = fraction2[good_fc]
        if test == "t":
            pvals = [ttest_ind(curr_data.iloc[cells1,i], curr_data.iloc[cells2,i])[1] for i in range(curr_data.shape[1])]
        elif test == "bimod":
            pvals = [differential_lrt(curr_data.iloc[cells1,i].values, curr_data.iloc[cells2,i].values) for i in range(curr_data.shape[1])]
        elif test == "wilcox":
            pvals = [ranksums(curr_data.iloc[cells1,i].values, curr_data.iloc[cells2,i].values)[1] for i in range(curr_data.shape[1])]
        d =  pd.DataFrame(data=np.array((pvals, log_fc, fraction1, fraction2)).T, columns=["pval", "log_fc", "pct.1", "pct.2"], index=curr_data.columns)
        return d 

    def compare_expression_between_groups(self, test="bimod", use_genes=None, use_raw=False):
        """
        Get log FC and P-value for each gene for each cluster between groups
        :return: dataframe containing each cluster and 
        """
        group_vals = np.unique(self._meta["group"].values)
        cluster_df = []
        ncells = []
        for c in np.unique(self._clusts):
            cells = self.get_cells_by_cluster(c, use_raw=use_raw)
            meta = self.get_metadata_by_cluster(c)
            cells0 = cells.iloc[get_subset_index(meta["group"], group_vals[0]),:]
            cells1 = cells.iloc[get_subset_index(meta["group"], group_vals[1]),:]
            if use_genes is not None:
                cells0 = cells0.loc[:,use_genes]
                cells1 = cells1.loc[:,use_genes]
            # log fold change for each gene for this cluster
            if use_raw:
                log_fc = np.array([np.log2((0.12+cells0.iloc[:,i].values.mean()))-np.log2(0.12+cells1.iloc[:,i].values.mean()) for i in range(cells0.shape[1])])
            else:
                log_fc = np.array([np.log1p(np.expm1(cells0.iloc[:,i].values).mean()) - np.log1p(np.expm1(cells1.iloc[:,i].values).mean()) for i in range(cells0.shape[1])])
            if test == "bimod":
                pvals = [differential_lrt(cells0.iloc[:,i].values, cells1.iloc[:,i].values) for i in range(cells1.shape[1])]
            elif test == "t":
                pvals = [ttest_ind(cells0.iloc[:,i].values, cells1.iloc[:,i].values)[1] for i in range(cells1.shape[1])]
            elif test == "wilcox":
                pvals = [ranksums(cells0.iloc[:,i].values, cells1.iloc[:,i].values)[1] for i in range(cells1.shape[1])]
            d = pd.DataFrame(data=np.array((cells0.mean(), cells1.mean(), pvals, log_fc, np.repeat(c, cells0.shape[1]).astype(np.int))).T, columns=["mean0", "mean1", "pval", "log_fc", "cluster"], index=cells0.columns)
            _,d["pval"],_,_ = multipletests(d["pval"], method="fdr_bh")
            d = d.sort_values(["pval", "log_fc"])
            cluster_df.append(d)
            ncells.append((cells0.shape[0], cells1.shape[0]))
        return pd.concat(cluster_df), ncells
    
    def plot_expression_between_groups(self,gene_names, test="bimod", plot_type="bar", figsize=(10,10), vmin=0, vmax=None, use_raw=False):
        """
        Plots expression between groups for sepecifc genes for all clusters
        """
        if not isinstance(gene_names, list):
            gene_names = [gene_names]

        group_vals = np.unique(self._meta["group"].values)
        cluster_df = []
        ncells = []
        f,ax = plt.subplots(nrows=len(gene_names), ncols=len(np.unique(self._clusts)), figsize=figsize)
        f.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)

        for i,g in enumerate(gene_names):
            for j,c in enumerate(np.unique(self._clusts)):
                cells = self.get_cells_by_cluster(c,use_raw=use_raw)
                # normalize to TPM
                meta = self.get_metadata_by_cluster(c)
                cells0 = cells.iloc[get_subset_index(meta["group"], group_vals[0]),:]
                cells1 = cells.iloc[get_subset_index(meta["group"], group_vals[1]),:]
                n0 = cells0.shape[0]
                n1 = cells1.shape[0]
                expr =np.hstack((cells0[g].values, cells1[g].values))

                ids = np.hstack((np.zeros((n0,)), np.ones((n1,))))
                temp = np.zeros_like(ids)
                d = pd.DataFrame(data=np.vstack((expr,ids)).T, columns=["expr", "group"])
                if len(gene_names) == 1:
                    curr_ax = ax[j]
                else:
                    curr_ax = ax[i][j]
                if plot_type is "bar":
                    sns.barplot(x="group", y="expr", data=d, ax=curr_ax, capsize=.2, errwidth=1) 
                elif plot_type is "violin":
                    sns.violinplot(x="group", y="expr", data=d, ax=curr_ax, capsize=.2, errwidth=1, palette="Set2_r",inner=None,linewidth=0)
                    sns.swarmplot(x="group", y="expr", data=d, ax=curr_ax, size=4, color=".3", linewidth=0)
                if test is "bimod":
                    pval = differential_lrt(cells0[g].values, cells1[g].values)
                if test is "wilcox":
                    pval = ranksums(cells0[g].values, cells1[g].values)[1]
                if test is "t":
                    pval = ttest_ind(cells0[g].values, cells1[g].values)[1]
                if vmax is not None:
                    curr_ax.set_ylim([vmin, vmax])
                curr_ax.set_title("Gene %s\nCluster %d\nP=%2.4f" % (g, c, pval))
                curr_ax.get_xaxis().set_visible(False)
                sns.despine(fig=f, ax=curr_ax, bottom=True, left=True)
                #sns.violinplo

    def volcano_plot(self, log_pval_thresh=5, log_fc_thresh=0.5, test="bimod", use_genes=None, use_raw=False):
        comparisons,ncells = self.compare_expression_between_groups(test=test, use_genes=use_genes, use_raw=use_raw)
        comparisons["pval"] = -np.log10(comparisons["pval"])
        ymax = comparisons["pval"].replace([np.inf, -np.inf], np.nan).dropna(how="all").max()
        n_clusts = len(np.unique(self._clusts))
        for i,c in enumerate(np.unique(self._clusts)):
            ax = plt.subplot(1,n_clusts,i+1)
            curr_vals = comparisons[comparisons["cluster"]==c]
            m = 6
            plt.plot(curr_vals["log_fc"], curr_vals["pval"], 'ko',markersize=m,markeredgewidth=0,linewidth=0)
            good_genes = curr_vals.loc[curr_vals["pval"]>log_pval_thresh,:]
            for g in good_genes.index:
                if g == "Egr2":
                    print("Clu=%d,Name=%s,logFC=%f,Pval=%f"%
                        (c,str(g),good_genes.loc[g,"log_fc"],good_genes.loc[g,"pval"]))
            for g in good_genes.index:
                x = good_genes.loc[g,"log_fc"]
                y = good_genes.loc[g,"pval"]
                plt.plot(x,y,'go',markersize=m,markeredgewidth=0,linewidth=0)
            good_genes = good_genes.loc[good_genes["log_fc"].abs()>log_fc_thresh,:]
            
            for g in good_genes.index:
                x = good_genes.loc[g,"log_fc"]
                y = good_genes.loc[g,"pval"]
                plt.plot(x,y,'ro',markersize=m,markeredgewidth=0,linewidth=0)
                plt.text(x,y,str(g),fontsize=18)
            plt.xlim([-2, 2])
            plt.ylim([-0.5, 1.2*ymax])
            ax.set_xticks([-2,0,2])
            if i > 0:
                ax.get_yaxis().set_visible(False)
            sns.despine()
            plt.tick_params(axis='both', which='major', labelsize=18)

    def dotplot_expression_across_clusters(self, gene_names,scale_max=500, cmap=plt.cm.viridis, clust_order=None):
        n_genes = len(gene_names)
        n_clusts = len(np.unique(self._clusts))
        uniq_clusts, clust_counts = np.unique(self._clusts, return_counts=True)
        if not clust_order:
            pos = range(n_clusts)
        avg = []
        num = []
        for i in range(n_genes):
            expr = self.get_expr_for_gene(gene_names[i], scaled=True).values        
            d = pd.DataFrame(np.array([expr, self._clusts]).T, columns=["expr", "cluster"])
            avg_expr = d.groupby("cluster").mean()
            avg_expr /= avg_expr.sum()
            avg_expr = avg_expr.values.flatten()
            num_expr = d.groupby("cluster").apply(lambda x: (x["expr"]>0).sum()).values.astype(np.float)
            num_expr /= clust_counts
            avg.append(avg_expr)
            num.append(num_expr)
        avg = np.vstack(avg)
        num = np.vstack(num)
        if clust_order:
            pos = []
            for i in range(n_genes):
                idx = np.argsort(-avg[i,:])
                for k in idx:
                    if k not in pos:
                        pos.append(k)
                        break
            print(avg.shape)
            print(pos)
            num = num[:,pos]
            avg = avg[:,pos]
            pos = range(num.shape[1])
        for i in range(n_genes):
            plt.scatter(pos, -i*np.ones_like(pos), s=num[i,:]*scale_max, c=avg[i,:], cmap=cmap, vmin=0, vmax=avg[i,:].max(),lw=0)
        plt.yticks(-np.array(range(n_genes)), gene_names)
        plt.axes().set_xticks(pos)
        plt.axes().set_xticklabels(pos)
        plt.xlabel('Cluster')

    def tsne_plot_orig_ident(self):
        plt.scatter(self._tsne[:,0], self._tsne[:,1], c=self._meta['orig_ident'], s=10,cmap=plt.cm.gist_rainbow,lw=0)

    def plot_expression_across_clusters(self, gene_names, plot_type="bar",figsize=None, clust_order=False,show_frame=True):
        """
        Plot array of genes x clusters
        """
        # so
        if clust_order is not None:
            clusts = self._clusts.copy()
            for i,c in enumerate(clust_order):
                clusts[self._clusts==c] = i
        else:
            clusts = self._clusts
        n_genes = len(gene_names)
        if figsize is None:
            f, ax = plt.subplots(n_genes, 1)
        else:
            f, ax = plt.subplots(n_genes, 1, figsize=figsize)
        f.tight_layout()
        for i in range(n_genes):
            expr = self.get_expr_for_gene(gene_names[i], scaled=False).values        
            #plt.title(name)
            d = pd.DataFrame(np.array([expr, clusts]).T, columns=["expr", "cluster"])
            if plot_type == "bar":
                sns.barplot(x="cluster", y="expr", data=d, ax=ax[i], capsize=.2, errwidth=1, color='gray') 
            elif plot_type == "violin":
                sns.violinplot(x="cluster", y="expr", data=d, ax=ax[i], capsize=.2, errwidth=1) 
            # get rid of the frame
            if not show_frame:
                for spine in ax[i].spines.values():
                    spine.set_visible(False)
                ax[i].get_xaxis().set_visible(False)
                ax[i].get_yaxis().set_visible(True)
                ax[i].set_ylabel(gene_names[i])
                if i == n_genes-1:
                    ax[i].tick_params(top='off', bottom='on', left='off', right='off', labelleft='off', labelbottom='on')
                else:
                    ax[i].tick_params(top='off', bottom='off', left='off', right='off', labelleft='off', labelbottom='on')


    def plot_bar_gene_expression(self, gene_names, nrow=None, ncol=None, ylim=None, figsize=(5,5),cmap=None):
        def _bar_plot(name,ylim=None,ax=None):
            expr = self.get_expr_for_gene(name, scaled=False).values        
            #plt.title(name)
            d = pd.DataFrame(np.array([expr, self._clusts]).T, columns=["expr", "cluster"])
            if cmap is None:
                sns.barplot(x="cluster", y="expr", data=d, ax=ax) 
            else:
                sns.barplot(x="cluster", y="expr", data=d, ax=ax,palette=cmap)
            if ax is None:
                ax = plt.axes()
            ax.set_title(name)
            if ylim is not None:
                ax.set_ylim([-1, ylim])
            sns.despine(ax=ax)
            ax.set_xlabel("")
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
        if isinstance(gene_names, list):
            if nrow is None and ncol is None:
                nrow = np.ceil(np.sqrt(len(gene_names)))
                ncol = nrow
            f,ax = plt.subplots(int(nrow), int(ncol),figsize=figsize)
            f.tight_layout()
            ax = np.array(ax).flatten()
            for i, name in enumerate(gene_names):
                _bar_plot(name,ylim,ax=ax[i])
        else:
            _bar_plot(gene_names,ylim)


    def plot_heatmap(self, gene_names,cmap=plt.cm.bwr,fontsize=16,use_imshow=False, ax=None, show_vlines=True):
        from scipy.stats.mstats import zscore
        if ax is None:
            ax = plt.axes()
        clust_sizes = [sum(self._clusts==i) for i in np.unique(self._clusts)]
        data = np.vstack([self._data.iloc[self._clusts==i,:].loc[:,gene_names].values for i in np.unique(self._clusts)]).T
        if use_imshow:
            ax.imshow(np.flipud(zscore(data,axis=1)), vmin=-2.5, vmax=2.5, cmap=cmap,interpolation='none', aspect='auto')
        else:
            ax.pcolor(np.flipud(zscore(data,axis=1)), vmin=-2.5, vmax=2.5, cmap=cmap)
        #plt.imshow(np.flipud(zscore(data,axis=1)),vmin=-2.5, vmax=2.5, cmap=cmap, aspect='auto', interpolation='none')
        ax.set_xlim([0,data.shape[1]])
        ax.set_ylim([0,data.shape[0]])
        if show_vlines:
            for i in np.cumsum(clust_sizes[:-1]):
                ax.axvline(i,color='k', linestyle='-')
        ax.set_yticks(np.arange(data.shape[0])+0.5)
        ax.set_yticklabels(gene_names[::-1],fontsize=fontsize)
        #ax.get_xaxis().set_fontsize(fontsize)
        plt.tick_params(axis='both', which='major', labelsize=fontsize)

    def plot_violin_gene_expression(self, gene_names, nrow=None, ncol=None, ylim=None):
        def _violin_plot(name,ylim=None):
            expr = self.get_expr_for_gene(name, scaled=False).values        
            plt.title(name)
            d = pd.DataFrame(np.array([expr, self._clusts]).T, columns=["expr", "cluster"])
            sns.violinplot(x="cluster", y="expr", data=d)
            if ylim is not None:
                plt.ylim([-1, ylim])

        if isinstance(gene_names, list):
            if nrow is None and ncol is None:
                nrow = np.ceil(np.sqrt(len(gene_names)))
                ncol = nrow
            for i, name in enumerate(gene_names):
                plt.subplot(nrow, ncol, i+1)
                _violin_plot(name,ylim)
        else:
            _violin_plot(gene_names,ylim)

    def plot_tsne_gene_expression(self, gene_names, scaled=True, nrow=None, ncol=None, s=10):
        """
        Plot the expression of a single gene in tSNE space
        """
        if isinstance(gene_names, list):            
            if nrow is None and ncol is None:
                nrow = np.ceil(np.sqrt(len(gene_names)))
                ncol = nrow
            for i,name in enumerate(gene_names):
                plt.subplot(nrow, ncol, i+1)
                expr = self.get_expr_for_gene(name, scaled=scaled)       
                plt.title(name,fontsize=16)
                plt.scatter(self._tsne[:,0], self._tsne[:,1],c=expr,cmap=plt.cm.jet,s=s,lw=0)
                plt.axis('off')
        else:
            expr = self.get_expr_for_gene(gene_names, scaled=scaled)        
            plt.title(gene_names,fontsize=16)
            plt.scatter(self._tsne[:,0], self._tsne[:,1],c=expr,cmap=plt.cm.jet,s=s,lw=0)
        #plt.axis('off')

def get_subset_index(p, i):
    return np.argwhere(p.values==i).flatten()

def differential_lrt(x,y,xmin=0):
    from scipy.stats import chi2
    lrtX = bimod_likelihood(x)
    lrtY = bimod_likelihood(y)
    lrtZ = bimod_likelihood(np.concatenate((x,y)))
    lrt_diff = 2 * (lrtX + lrtY - lrtZ)

    return chi2.pdf(x=lrt_diff, df=3)

def min_max(x, x_min, x_max):
    x[x>x_max] = x_max
    x[x<x_min] = x_min
    return x

def bimod_likelihood(x, xmin=0):
    x1 = x[x <= xmin]
    x2 = x[x > xmin]
    #xal = MinMax(x2.shape[0]/x.shape[0], 1e-5, 1-1e-5)
    xal = len(x2)/float(len(x))
    if xal < 1e-5:
        xal = 1e-5
    if xal > 1-1e-5:
        xal = 1-1e-5
    likA = x1.shape[0] * np.log(1-xal)
    if len(x2) < 2:
        mysd = 1
    else:
        mysd = np.std(x2)
    likB = x2.shape[0]*np.log(xal) + np.sum(np.log(norm.pdf(x2, np.mean(x2), mysd)))
    return likA + likB

def top_genes_per_clust(markers, N=3, pval_thresh=1e-6, return_unique=False):
    top_genes = []
    clusts = np.unique(markers["cluster"])
    for c in clusts:
        curr_genes = markers[markers["cluster"] == c]
        curr_genes = curr_genes[curr_genes["pval"] < pval_thresh]
        top_genes.extend(list(curr_genes.index[:N]))
    if return_unique:
        dup_top_genes = top_genes
        top_genes = []
        for i in dup_top_genes:
            if i not in top_genes:
                top_genes.append(i)
        #top_genes = list(set(top_genes))
    return top_genes
