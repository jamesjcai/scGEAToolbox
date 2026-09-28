function [c] = ml_SnnDpc(s, cluK, knnK)
% Clustering cell embeddings using SNNDPC - a SNN clustering algorithm
%
% http://mlwiki.org/index.php/SNN_Clustering#SSN_Clustering_Algorithm
% https://link.springer.com/article/10.1007/s12539-019-00357-4
%
% The work is done by PKG.E_SNNDPC, which reproduces the reference
% implementation in EXTERNAL/ML_SNNDPC/SNNDPC_ORI exactly but in O(N*K)
% rather than O(N^2) time and memory. SNNDPC_ORI is kept only as the
% reference TESTS/SNNDPCLOUVAINTEST checks against.

if nargin < 3, knnK = 4; end
if nargin < 2, cluK = 10; end

c = pkg.e_snndpc(s, cluK, knnK);

end
