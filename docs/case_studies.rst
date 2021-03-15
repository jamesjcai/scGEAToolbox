.. _case_studies:

Case Studies and Tutorials
==========================

Download data files from GEO
----------------------------

From GEO database, we obtain the FTP links to the data files we need. Here we use a data set from sample GSM3535276 as an example ( https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3535276). The sample is human AXLN1 lymphatic endothelial cells.

.. raw:: html

  <table cellpadding="2" cellspacing="2" width="600"><tr bgcolor="#eeeeee" valign="top"><td align="middle" bgcolor="#CCCCCC"><strong>Supplementary file</strong></td>
  <td align="middle" bgcolor="#CCCCCC"><strong>Size</strong></td>
  <td align="middle" bgcolor="#CCCCCC"><strong>Download</strong></td>
  <td align="middle" bgcolor="#CCCCCC"><strong>File type/resource</strong></td>
  </tr>
  <tr valign="top"><td bgcolor="#DEEBDC">GSM3535276_AXLN1_barcodes.tsv.gz</td>
  <td bgcolor="#DEEBDC" title="34376">33.6 Kb</td>
  <td bgcolor="#DEEBDC"><a href="ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3535nnn/GSM3535276/suppl/GSM3535276%5FAXLN1%5Fbarcodes%2Etsv%2Egz">(ftp)</a><a href="/geo/download/?acc=GSM3535276&amp;format=file&amp;file=GSM3535276%5FAXLN1%5Fbarcodes%2Etsv%2Egz">(http)</a></td>
  <td bgcolor="#DEEBDC">TSV</td>
  </tr>
  <tr valign="top"><td bgcolor="#EEEEEE">GSM3535276_AXLN1_genes.tsv.gz</td>
  <td bgcolor="#EEEEEE" title="257278">251.2 Kb</td>
  <td bgcolor="#EEEEEE"><a href="ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3535nnn/GSM3535276/suppl/GSM3535276%5FAXLN1%5Fgenes%2Etsv%2Egz">(ftp)</a><a href="/geo/download/?acc=GSM3535276&amp;format=file&amp;file=GSM3535276%5FAXLN1%5Fgenes%2Etsv%2Egz">(http)</a></td>
  <td bgcolor="#EEEEEE">TSV</td>
  </tr>
  <tr valign="top"><td bgcolor="#DEEBDC">GSM3535276_AXLN1_matrix.mtx.gz</td>
  <td bgcolor="#DEEBDC" title="47972933">45.8 Mb</td>
  <td bgcolor="#DEEBDC"><a href="ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3535nnn/GSM3535276/suppl/GSM3535276%5FAXLN1%5Fmatrix%2Emtx%2Egz">(ftp)</a><a href="/geo/download/?acc=GSM3535276&amp;format=file&amp;file=GSM3535276%5FAXLN1%5Fmatrix%2Emtx%2Egz">(http)</a></td>
  <td bgcolor="#DEEBDC">MTX</td>
  </tr>
  </table>


We can use `gunzip` function directly download and unzip the files.

.. code-block:: matlab

  gunzip('https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3535nnn/GSM3535276/suppl/GSM3535276_AXLN1_barcodes.tsv.gz')
  gunzip('https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3535nnn/GSM3535276/suppl/GSM3535276_AXLN1_genes.tsv.gz');
  gunzip('https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3535nnn/GSM3535276/suppl/GSM3535276_AXLN1_matrix.mtx.gz');

We can then use the code below to import data into `MATLAB`.

.. code-block:: matlab

  [X,genelist,barcodelist]=sc_readmtxfile('matrix.mtx','features.tsv','barcodes.tsv',2);
  [X,genelist]=sc_qcfilter(X,genelist);
  [X,genelist]=sc_selectg(X,genelist,1,0.05);
  sc_scatter(X,genelist)



Process 10x Genomics raw data
-----------------------------
In the 10x Genomics folder, there are three files, namely, matrix.mtx, features.tsv and barcodes.tsv. Here is the best practice of raw data processing.

.. code-block::
  
  [X,genelist,barcodelist]=sc_readmtxfile('matrix.mtx','features.tsv','barcodes.tsv',2);
  [X,genelist]=sc_qcfilter(X,genelist);
  [X,genelist]=sc_selectg(X,genelist,1,0.05);
  [Xnorm]=sc_norm(X,'type','deseq');
  [~,Xhvg]=sc_hvg(Xnorm,genelist,true);
  [s_tsne]=sc_tsne(Xhvg(1:2000,:),3,false,false);
  sc_scatter(X,genelist,s_tsne)

Overview of `SC_SCATTER`
-----------------------
`SC_SCATTER` can be used to visualize `SCE` class/object. Below are links to several case studies and examples using sc_scatter_sce to explore high-dimensional data. All examples are below are publically available through GitHub.

|Overview of sc_scatter|

Import Seurat RData
-------------------
For example, we are trying to read files from `https://www.synapse.org/#!Synapse:syn22855256 <https://www.synapse.org/#!Synapse:syn22855256>`_. They are described as `pbmc_discovery_v1.RData` and `pbmc_replication_v1.RData` are Seurat objects containing the gene expression raw counts and log normalized data, the phenotype Label ("CI" for MCI, "C" for control) and the inferred cell identity of the discovery and replication cohort, respectively. 

.. code-block:: r

  library(Seurat)
  library(Matrix)
  load('pbmc_discovery_v1.RData')
  countMatrix <- pbmc_discovery@assays$RNA@counts
  writeMM(obj = countMatrix, file = 'matrix.mtx')
  writeLines(text = rownames(countMatrix), con = 'features.tsv')
  writeLines(text = colnames(countMatrix), con = 'barcodes.tsv')
  metadata <- pbmc_discovery@meta.data
  write.csv(x = metadata, file = 'metadata.csv', quote = FALSE)

After exporting Seurate object data into the three files, you can then use MATLAB to read the files:
  
.. code-block:: matlab

  [X,genelist,barcodelist]=sc_readmtxfile('matrix.mtx','features.tsv','barcodes.tsv',1);  
  sce=SingleCellExperiment(X,genelist);
  T=readtable('metadata.csv')
  c=string(T.Label);
  sce.c_batch_id=c;
  sc_scatter(sce)
  
.. |Overview of sc_scatter| image:: https://github.com/jamesjcai/scGEAToolbox/raw/master/resources/Tooltips.png
   :target: https://github.com/jamesjcai/scGEAToolbox/raw/master/resources/Tooltips.png
  
