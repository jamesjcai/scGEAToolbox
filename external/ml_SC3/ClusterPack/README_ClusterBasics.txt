----------------------------------------------------------------------
CLUSTERBASICS   README    Alexander Strehl    Version 2.0   2011-05-01
----------------------------------------------------------------------
This package contains MATLAB implementation of a variety of cluster
functions as described in 
      A. Strehl, J. Ghosh and R. Mooney, "Impact of Similarity
      Measures on Web-page Clustering", Proc. of the 17th National
      Conference on Artificial Intelligence: Workshop of Artificial
      Intelligence for Web Search (AAAI 2000), July 2000, Austin,
      Texas, pp. 58-64 
and
      A. Strehl and J. Ghosh, "Relationship-based Clustering and
      Visualization for High-dimensional Data Mining", INFORMS
      Journal on Computing, pages 208-230, Spring 2003
----------------------------------------------------------------------
To use it copy all 31 files
      ClusterBasicsTest.m
      README_ClusterBasics.txt
      checkcl.m
      checks.m
      clagglmin.m
      clcgraph.m
      clgraph.m
      clhgraph.m
      clkmeans.m
      clrand.m
      clucent.m
      cmetis.m
      evalbalance.m
      evalf.m
      evalmse.m
      evalmutual.m
      fastchangem.m
      hmetis.m
      ints.m
      kmeans.m
      mapasc.m
      mapdense.m
      metis.m
      onetomax.m
      sgraph.m
      simbjac.m
      simcorr.m
      simcosi.m
      simeucl.m
      simxjac.m
      wgraph.m
to a directory in your MATLAB path or the current directory in
your MATLAB environment. Also, make sure that
      pmetis
      shmetis
or
      pmetis.exe
      shmetis.exe
are in your machine's path and executable. Metis is a separate
program available from
      http://www-users.cs.umn/edu/~karypis/metis/download.html
The Linux & Win32 binaries are included with this distribution. Then,
in MATLAB just type
      ClusterBasicsTest;
to run a test of all the functions. Type
      help ClusterBasicsTest 
to obain other info or
      edit ClusterBasicsTest
to see the function calls performed for the test.
All functions have been tested on WIN7 Octave 3.2.4 and LNX86 Matlab 5.2.0.3084.
----------------------------------------------------------------------
The following functions are intended for user invocation. They are
arranged into four groups of activities:
* Similarity Matrix Computation. The following functions compute
a matrix containing all pairwise similarities between the row
vectors in a and b:
      simeucl(a,b)
      simcosi(a,b)
      simcorr(a,b)
      simbjac(a,b)
      simxjac(a,b)
* Clustering Algorithms. The following functions compute the
cluster labels from 1 to k for all the objects described by row
vectors in x using the similarity semantic sfct (e.g., 'simeucl'):
      clagglmin(x,k,sfct)
      clcgraph(x,k,sfct)
      clgraph(x,k,sfct)
      clhgraph(x,k,sfct)
      clkmeans(x,k,sfct,oldcl)
      clrand(x,k,sfct)
* Cluster Quality Evaluation Functions. The following functions
evaluate the quality of a clustering cl given a human imposed
categorization trueclass for the data matrix x using similarity
semantic sfct:
      evalbalance(trueclass,cl,x,sfct)
      evalf(trueclass,cl,x,sfct)
      evalmse(trueclass,cl,x,sfct)
      evalmutual(trueclass,cl,x,sfct)
* Data Cleaning. The following functions can be used to
check and correct problems with a clustering cl or a similarity
matrix s:
      checkcl(cl)
      checks(s)
Please note that not all functions use all parameters.  Please type
'help functionname' for more information on an individual function.
----------------------------------------------------------------------
For questions, comments & services please contact
      alexander@strehl.com
This package and potential updates are available from
      http://strehl.com
----------------------------------------------------------------------
Finally, copyright (c) 1998-2011 by Alexander Strehl. All rights
reserved.  License is granted to copy, to use, and to make and to use
derivative works for research purposes, provided that the Alexander
Strehl copyright notice and this license notice is included in all
copies and any derivatives works and in all related
documentation. Alexander Strehl grants no other licenses expressed or
implied and the licensee acknowleges that Alexander Strehl has no
liability for licensee's use or for any derivative works made by
licensee. This software is provided as is. Alexander Strehl disclaims
and licensee agrees that all warranties, express or implied, including
without limitation the implied warranties of merchantability and
fitness for a particular purpose. Notwithstanding any other provision
contained herein, any liability for damages resulting from the
software or its use is expressly disclaimed, including consequential
or any other indirect damages, whether arising in contract, tort
(including negligence) or strict liability, even if Alexander Strehl
is advised of the possibility of such damages. Enjoy!
----------------------------------------------------------------------


