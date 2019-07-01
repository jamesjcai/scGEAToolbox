Network Enhancement: a general method to denoise weighted biological networks
==============================================================================

Project website: [http://snap.stanford.edu/ne](http://snap.stanford.edu/ne).

Overview
--------

Networks are abundant in many areas of biology. These networks often entail non-trivial 
topological features and patterns critical to understanding interactions within the 
natural system. However, networks observed in real-world are typically noisy. The presence 
of high levels of noise can hamper discovery of structures and dynamics present in the network. 

We propose *Network Enhancement (NE),* a novel method for improving the 
signal-to-noise ratio of a symmetric networks and thereby facilitating the downstream network 
analysis. NE leverages the transitive edges of a network by exploiting local structures to 
strengthen the signal within clusters and weaken the signal between clusters. At the same time 
NE also alleviates the corrupted links in the network by imposing a normalization that removes 
weak edges by enforcing sparsity. NE is supported by theoretical justifications for its 
convergence and performance in improving community detection outcomes. 

The method provides theoretical guarantees as well as excellent empirical performance on 
many biological problems. The approach can be incorporated into any weighted network analysis 
pipeline and can lead to improved downstream analysis. 


Running NE in Matlab
--------------------

At current stage, we provide examples showing how to apply NE to two problems in biology, which 
are discussed in the manuscript. All datasets required to run the examples are included in this
repository. 

1) The first example is about fine-grained species identification. Run the example as:

```
    run_butterfly_network.m 
```

This script reports retrieval accuracy values and generates a retrieval curve for the task of 
butterfly species identification, as reported in the manuscript.

2) The second example is about denoising Hi-C interaction networks. We provide sample data for 
Chrom 16 with two measurement resolutions, 1kb and 5kb Hi-C data. Run the example as:

``` 
    run_hiC_network.m
```
    
In order to use community detection with the Hi-C interaction networks, you need to 
compile two C++ files. Instructions are provided in script ``run_hiC_network.m``.
