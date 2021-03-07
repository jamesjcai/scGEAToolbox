.. _case_studies:

Case Studies and Tutorials
==========================
`sc_scatter_sce` can be used to visualize `SCE` class/object. Below are links to several case studies and examples using sc_scatter_sce to explore high-dimensional data. All examples are below are publically available through GitHub.

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

|Overview of sc_scatter|

Using SC_SCATTER_SCE to explore scRNA-seq data stored as SCE class
------------------------------------------------------------------
.. raw:: html

         <div style="position: relative; padding-bottom: 10px; height: 0; overflow: hidden; max-width: 100%; height: auto;">

         <iframe width="560" height="315" src="https://www.youtube.com/embed/6IYl35dTwGw" frameborder="0" allow="autoplay; encrypted-media" allowfullscreen></iframe>
         </div>


Label cell type interactively with SC_SCATTER_SCE
-------------------------------------------------
.. raw:: html

         <div style="position: relative; padding-bottom: 10px; height: 0; overflow: hidden; max-width: 100%; height: auto;">

         <iframe width="560" height="315" src="https://www.youtube.com/embed/HRQiXX3Jwpg" frameborder="0" allow="autoplay; encrypted-media" allowfullscreen></iframe>
         </div>

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
  <tr><td><a href="/Traces/study/?acc=SRX5197633">SRA Run Selector</a></td></tr>
  <tr><td class="message">Raw data are available in SRA</td></tr>
  <tr><td class="message">Processed data provided as supplementary file</td></tr>
  </table>

We can use `gunzip` function directly download and unzip the files

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



.. |Overview of sc_scatter| image:: https://pbs.twimg.com/media/ErAPaabW4AEwAvf?format=png&name=large
   :target: https://youtu.be/6IYl35dTwGw
   
