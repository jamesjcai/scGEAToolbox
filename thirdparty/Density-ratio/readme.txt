These functions implement two approaches based on density-ratio for discovering 
clusters with varying densities, proposed by Ye Zhu, Kai Ming Ting, and Mark J. Carman. 
"Density-ratio based clustering for discovering clusters with varying densities." 
Pattern Recognition (2016). http://www.sciencedirect.com/science/article/pii/S0031320316301571

Written by Ye Zhu, Monash University, July 2016, version 1.1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Abstract of the paper

Density-based clustering algorithms are able to identify clusters of arbitrary shapes and sizes in a dataset
which contains noise. It is well-known that most of these algorithms, which use a global density
threshold, have difficulty identifying all clusters in a dataset having clusters of greatly varying densities.
This paper identifies and analyses the condition under which density-based clustering algorithms fail in
this scenario. It proposes a density-ratio based method to overcome this weakness, and reveals that it can
be implemented in two approaches. One approach is to modify a density-based clustering algorithm to
do density-ratio based clustering by using its density estimator to compute density-ratio. The other
approach involves rescaling the given dataset only. An existing density-based clustering algorithm, which
is applied to the rescaled dataset, can find all clusters with varying densities that would otherwise impossible 
had the same algorithm been applied to the unscaled dataset.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

There are three scripts you can run individually:

testReCon.m and testRescale.m demostrates the ReCon approach and ReScale approach, respectively. 

VisualiseReScale.m compares the MDS plots between on the original and ReScaled datasets.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

This code was used for the following publication (Bibtex format):

@article{zhu2016density,
  title={Density-ratio based clustering for discovering clusters with varying densities},
  author={Zhu, Ye and Ting, Kai Ming and Carman, Mark J},
  journal={Pattern Recognition},
  volume={60},
  pages={983--997},
  year={2016},
  publisher={Elsevier}
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

References: 
[1] M. Ankrest, M. Breunig, H. Kriegel, J. Sander, 
OPTICS: Ordering Points To Identify the Clustering Structure, 
available from www.dbs.informatik.uni-muenchen.de/cgi-bin/papers?query=--CO
[2] M. Daszykowski, B. Walczak, D.L. Massart, Looking for natural  
patterns in analytical data. Part 2. Tracing local density 
with OPTICS, J. Chem. Inf. Comput. Sci. 42 (2002) 500-507
[3] M. Ester, H. peter Kriegel, J. S, X. Xu, A density-based algorithm for 
discovering clusters in large spatial databases with noise, in: Proceedings 
of the Second International Conference on Knowledge Discovery and Data Mining 
(KDD-96), AAAI Press, 1996, pp. 226-231.
[4] L. Ertoz, M. Steinbach, V. Kumar, Finding clusters of different sizes, 
shapes, and densities in noisy, high dimensional data., in: SDM, 2003, pp. 47-58.
