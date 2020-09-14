% SCGEATOOLBOX
%
% Files
%   cdgea                     - - Changes to GEAToolbox directory
%   demo_script1              - Demonstration of Filter, Normalization and Batch Correction of Data in scGEAToolbox
%   demo_script2              - Demonstration of Feature Selection Functions in scGEAToolbox
%   demo_script3              - Demonstration of Visualization Functions in scGEAToolbox
%   demo_script4              - Demonstration of Clustering Functions in scGEAToolbox
%   demo_script5              - Demonstration of Pseudotime Analysis and Gene Network Functions in scGEAToolbox
%   demo_script6              - Demonstration of DE Analysis and Marker Gene Identification Functions in scGEAToolbox
%   fun_cmp_clusters          - https://academic.oup.com/bib/advance-article/doi/10.1093/bib/bby076/5077112
%   fun_num_cluster           - 
%   fun_sim_matrix            - compute cell-to-cell similarity matrix
%   i_get_hsmm_tf             - 
%   i_go_analysis             - http://geneontology.org/docs/go-annotation-file-gaf-format-2.0/
%   i_gscatter3               - 
%   i_labelcluster            - 
%   i_myscatter               - 
%   i_plot_pseudotimeseries   - 
%   i_pseudotime_by_splinefit - s_selected=s_phate0g1(cell_i,:);
%   i_quickload_examples      - 
%   i_stemscatter             - 
%   mminfo                    - function  [rows, cols, entries, rep, field, symmetry] = mminfo(filename)
%   mmread                    - function  [A] = mmread(filename)
%   mmwrite                   - Function: mmwrite(filename,A,comment,field,precision)
%   mnncorrect                - https://rdrr.io/bioc/scran/src/R/mnnCorrect.R
%   myboxplot                 - 
%   norm_deseq                - For DESeq normalization, the geometric mean for each gene was computed 
%   norm_libsize              - 
%   run_busseq                - NCELLS=2000; NGENES=400;
%   run_celltypeassignation   - A custom R code using fGSEA to assign cell type from given ranked list of
%   run_combat                - 
%   run_combat2               - 
%   run_combatn               - 
%   run_csndm                 - 
%   run_enrichr               - Run Enrichr
%   run_fgsea                 - 
%   run_fitsne                - RUN_PHATE 
%   run_g2s3                  - 
%   run_genie3                - RUN_GENIE3
%   run_glmpca                - 
%   run_gorilla               - Run GOrilla
%   run_gsea                  - 
%   run_hcp                   - Run HCP (Hidden Covariates with Prior) to normalize RNA-seq data
%   run_magic                 - 
%   run_mast                  - 
%   run_mcimpute              - 
%   run_monocle               - 
%   run_pcnet                 - if nargin<2, plotit=false; end
%   run_phate                 - 
%   run_scode                 - if nargin<2, plotit=false; end
%   run_simlr                 - - 
%   run_singler               - 
%   run_sinnlrr               - SinNLRR - 
%   run_soptsc                - run_soptsc - 
%   run_umap                  - 
%   sc_alona                  - sc_alona - 
%   sc_cdr                    - - calcualtes Cellular Detection Rate (CDR) of cells
%   sc_cellscatter            - 
%   sc_celltypecaller         - https://academic.oup.com/database/article/doi/10.1093/database/baz046/5427041
%   sc_celltypeexplorer       - load cleandata.mat
%   sc_celltypeexplorer_auto  - 
%   sc_cluster                - 
%   sc_clustshow              - if min(size(s))>3, error('S is coordinates of dimensional reduction.'); end
%   sc_deg                    - https://satijalab.org/seurat/v3.1/de_vignette.html
%   sc_diffuse                - Outputs:
%   sc_diptest                - USAGE:
%   sc_dotplot                - 
%   sc_embed4w                - 
%   sc_explorer               - 
%   sc_filterc                - 
%   sc_filterg                - 
%   sc_grnetwork              - 
%   sc_heg                    - HEGs - highly expressed genes
%   sc_hvg                    - HVGs selection - This method uses the CV^2 on normalized count data to 
%   sc_lcod                   - Lowess coefficient of dispersion (LCOD) analysis
%   sc_marker                 - 
%   sc_markerexplorer         - load cleandata.mat
%   sc_markerscatter          - SC_MARKERSCATTER(X,genelist,g,s,methodid)
%   sc_merge2data             - 
%   sc_mmread                 - 
%   sc_norm                   - 
%   sc_pcnet                  - [A]=sc_pcnet(X,ncom)      % X = expression matrix of genes x cells
%   sc_pcnetpar               - [A]=sc_pcnetpar(X,ncom)
%   sc_pickmarkers            - IDV - cluster ids of cells
%   sc_pickmarkers2           - 
%   sc_pseudotimeexplorer     - load cleandata.mat
%   sc_qcfilter               - 
%   sc_qcfilter2              - 
%   sc_qcmetrics              - {
%   sc_readfile               - 
%   sc_readh5file             - https://www.mathworks.com/help/matlab/hdf5-files.html
%   sc_readmtxfile            - 
%   sc_readtsvfile            - 
%   sc_rmmtcells              - 
%   sc_rmmtgenes              - 
%   sc_sc3                    - SC3 - consensus clustering of single-cell RNA-seq data
%   sc_scatter                - 
%   sc_scatter3               - 
%   sc_scatter3x              - 
%   sc_scattermarker          - SC_SCATTERMARKER(X,genelist,g,s,methodid)
%   sc_sct                    - https://satijalab.org/seurat/v3.0/sctransform_vignette.html
%   sc_selectc                - 
%   sc_selectg                - https://www.nature.com/articles/s41592-018-0254-1
%   sc_similarg               - 
%   sc_splinefit              - - 
%   sc_splinefit2             - 
%   sc_stat                   - 
%   sc_stat2                  - 
%   sc_stem3                  - 
%   sc_stemscatter            - 
%   sc_trajectory             - 
%   sc_transform              - 
%   sc_tscan                  - ref: https://academic.oup.com/nar/article/44/13/e117/2457590
%   sc_tsne                   - see also: RUN_PHATE, RUN_UMAP
%   sc_veg                    - 
%   sc_view2grpcells          - 
%   sc_viewngrpcells          - 
%   sc_writefile              - 
%   splinefit                 - Fit a spline to noisy data.
%   zip_redistribution        - 


