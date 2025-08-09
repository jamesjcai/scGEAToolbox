.. _case_studies:

Case Studies and Tutorials
==========================

Download 10x Genomics data files from GEO
-----------------------------------------

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

  gunzip('https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3535nnn/GSM3535276/suppl/GSM3535276_AXLN1_matrix.mtx.gz');
  gunzip('https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3535nnn/GSM3535276/suppl/GSM3535276_AXLN1_genes.tsv.gz');

We can then use the code below to import data into `MATLAB`.

.. code-block:: matlab

  [X,g]=sc_readmtxfile('GSM3535276_AXLN1_matrix.mtx','GSM3535276_AXLN1_genes.tsv');
  sce=SingleCellExperiment(X,g);
  scgeatool(sce)


Process downloaded 10x Genomics data files
------------------------------------------
In a 10x Genomics data folder, there should be matrix.mtx and genes.tsv. Here is the commandline code for raw data processing.

.. code-block::
  
  [X,g]=sc_readmtxfile('matrix.mtx','genes.tsv');
  [X,g]=sc_qcfilter(X,g);
  [X,g]=sc_selectg(X,g,1,0.05);
  [s]=sc_tsne(X);
  sce=SingleCellExperiment(X,g,s);
  scgeatool(sce)


Download Drop-seq data files from GEO
-------------------------------------

From GEO database, we obtain the FTP links to the data files we need. Here we use a data set from sample GSM3036814 as an example (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3036814). The sample is mouse lung cells.

.. raw:: html

  <table cellpadding="2" cellspacing="2" width="600"><tr bgcolor="#eeeeee" valign="top"><td align="middle" bgcolor="#CCCCCC"><strong>Supplementary file</strong></td>
  <td align="middle" bgcolor="#CCCCCC"><strong>Size</strong></td>
  <td align="middle" bgcolor="#CCCCCC"><strong>Download</strong></td>
  <td align="middle" bgcolor="#CCCCCC"><strong>File type/resource</strong></td>
  </tr>
  <tr valign="top"><td bgcolor="#DEEBDC">GSM3036814_Control_6_Mouse_lung_digital_gene_expression_6000.dge.txt.gz</td>
  <td bgcolor="#DEEBDC" title="1797463">1.7 Mb</td>
  <td bgcolor="#DEEBDC"><a href="ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3036nnn/GSM3036814/suppl/GSM3036814%5FControl%5F6%5FMouse%5Flung%5Fdigital%5Fgene%5Fexpression%5F6000%2Edge%2Etxt%2Egz">(ftp)</a><a href="/geo/download/?acc=GSM3036814&amp;format=file&amp;file=GSM3036814%5FControl%5F6%5FMouse%5Flung%5Fdigital%5Fgene%5Fexpression%5F6000%2Edge%2Etxt%2Egz">(http)</a></td>
  <td bgcolor="#DEEBDC">TXT</td>
  </tr>
  </table>
  
We can use `gunzip` function directly download and unzip the files.

.. code-block:: matlab

  gunzip('https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3036nnn/GSM3036814/suppl/GSM3036814_Control_6_Mouse_lung_digital_gene_expression_6000.dge.txt.gz')
  
We can then use the code below to import data into `MATLAB`.

.. code-block:: matlab

  [X,g]=sc_readtsvfile('GSM3036814_Control_6_Mouse_lung_digital_gene_expression_6000.dge.txt');
  [X,g]=sc_qcfilter(X,g);
  [X,g]=sc_selectg(X,g,1,0.05);
  [s]=sc_tsne(X);
  sce=SingleCellExperiment(X,g,s);
  scgeatool(sce)

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
  scgeatool(sce)

Import data from a TSV/Excel file
---------------------------------
If your scRNA-seq data is in Excel file, save it as TSV or CSV a file with the format like this:

.. code-block:: text

  genes	X1	X2	X3	X4	X5	X6	X7	X8	X9
  NOC2L	1	1	2	3	3	2	0	1	3	
  HES4	50	15	19	50	8	87	23	25	29
  ISG15	279	312	425	180	406	408	335	403	398
  AGRN	3	4	9	5	2	3	8	8	9	
  SDF4	2	2	4	0	5	0	4	2	5	
  B3GALT6	2	1	0	0	1	0	1	1	0	
  UBE2J2	1	2	3	1	1	1	6	3	4	
  SCNN1D	0	1	0	0	0	0	0	0	0	
  ACAP3	1	3	1	0	1	0	0	1	0

Then you can use function `sc_readtsvfile` to import the data. Here is an example:

.. code-block:: matlab

  cdgea;
  [X,g]=sc_readtsvfile('example_data\GSM3204304_P_P_Expr.csv');

Visualize data in 6D
--------------------

.. code-block:: matlab

  cdgea;
  load example_data\example10xdata.mat
  % s=sc_tsne(X,6,false,true);
  s=s_tsne6;    % using pre-computed 6-d embedding S_TSNE6
  gui.sc_multiembeddings(s(:,1:3),s(:,4:6));
  
Here is what you should get:
  
|sixdview|

.. |sixdview| image:: https://github.com/jamesjcai/scGEAToolbox/raw/main/resources/Images/six_d.png
   :width: 250
   :target: https://github.com/jamesjcai/scGEAToolbox/raw/main/resources/Images/six_d.png
  
