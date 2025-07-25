# CSN
Cell-specific Network Constructed by Single-cell RNA Sequencing Data

    function ndm = csndm(data,alpha,boxsize,normalize)

 Construction of network degree matrix
 
 The function performs the transformation from gene expression matrix to network degree matrix (ndm).
 
 data: Gene expression matrix (TPM/FPKM/RPKM/count), rows = genes, columns = cells
 
 alpha: Significant level (eg. 0.001, 0.01, 0.05 ...), Default = 0.01
 
 boxsize: Size of neighborhood, Default = 0.1
 
 normalize: 1: result is normalized (Default); 0: result is not normalized
 
    
    
    
 
    function csn = csnet(data,c,alpha,boxsize,weighted)
 Construction of cell-specific network
 
 The function performs the transformation from gene expression matrix to cell-specific network (csn).
 
 data: Gene expression matrix, rows = genes, columns = cells
 
 c: Construct the CSNs for all cells, set c = [] (Default); Construct the CSN for cell k, set c = k
 
 alpha: Significant level (eg. 0.001, 0.01, 0.05 ...), larger alpha leads to more edges, Default = 0.01
 
 boxsize: Size of neighborhood, Default = 0.1
 
 weighted: 1: edge is weighted; 0: edge is not weighted (Default)
 
 csn: Cell-specific network, the kth CSN is in csn{k}, rows = genes, columns = genes
 
 Note that too many cells or genes may lead to out of memory.
 
    
    
    
 
    function edge = csnedge(gx,gy,boxsize)

 The normalized statistic of edge x-y
 
 gx gy: Gene expression values of gene x and gene y. If there are n cells, gx and gy are 1-by-n vectors.
 
 boxsize: Size of neighborhood, Default = 0.1
 
 edge: 1-by-n vector, the normalized statistic of edge x-y in all cells
