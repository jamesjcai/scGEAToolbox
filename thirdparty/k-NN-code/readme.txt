
Authors: Jose Costa, Alfred Hero

Date: April, 2004



0. Copyright notice

 Matab scripts for intrinsic dimension and entropy estimation using k-nearest
 neighbor graphs. The details of the algorithms can be found in:

   J. A. Costa and A. O Hero, "Entropic Graphs for Manifold Learning",
   Proc. of IEEE Asilomar Conf. on Signals, Systems, and Computers,
   Pacific Groove, CA, November, 2003.

   J. A. Costa and A. O. Hero, "Geodesic Entropic Graphs for Dimension and
   Entropy Estimation in Manifold Learning", 
   to appear in IEEE Trans. on Signal Processing, Aug., 2004. 

 Published reports of research using the code provided here (or a modified version)
 should cite the two articles referenced above.

 Comments and questions are welcome. We would also appreciate hearing about 
 how you used this code, improvements made to it, etc. You are free to modify the
 code, as long as you reference the original contributors.




1. List of Matlab files:

swiss_roll_example.m - example of using k-NN functions
                       for dimension and entropy estimation (on the swiss roll)


knn_graph_estim_1.m - function that estimates the intrinsic dimension and entropy
                      of a data set. Uses a fast C implementation (kNNlengthmex.c)
                      based on KD trees


knn_graph_estim_2.m - function that estimates the intrinsic dimension and entropy
                      of a data set. Uses an obvious (knn_distmat.m) implementation


The choice to use knn_graph_estim_1.m or knn_graph_estim_2.m comes from the extrinsic
dimensionality of the space (d) and the number of samples (n). The KD implementation of
the k-NN graph has complexity O(k d n log(n)), while the obvious implementation has
complexity O(k n^2 log(n)). For low dimensional spaces (i.e., d<<n), knn_graph_estim_1.m runs
much faster, but for high dimensional spaces (i.e., d>n), knn_graph_estim_2.m will be
faster. This is the case of spaces of (high-resolution) images.


plot_knn.m - (very unpolished) function to visualize 2D or 3D k-NN graphs


estim_beta_k_NN.m - simple script to estimate the constants necessary for entropy estimation




2. Before starting

Before running these scripts, you should compile the C libraries in the Matlab command
window as usual, e.g.:

mex kNNlengthmex.c
mex kNNgraphmex.c

This code was tested on Windows XP and Linux systems, using Matlab R12 and R13.



Jose Costa (jcosta@umich.edu)

