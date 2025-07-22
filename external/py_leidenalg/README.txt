Installation: 
conda env create -f environment.yml
conda activate scanpy_env_311 


ADD TO MATLAB PATH:  
HOME -> Set Path -> Add Folder (select this leiden_clustering_sparse directory)


Change python intepreter direction for matlab by:
In Terminal write 'conda env list' and select the directory shown for 'leiden_clustering'. 
EXAMPLE
(base) PS C:\Users\selim\Documents\vs_code_working_dir> conda env list
# conda environments:
#
base                     F:\Anaconda
scanpy_env_311        *  F:\Anaconda\envs\leiden_clustering
scanpy_env               F:\Anaconda\envs\scanpy_env


IMPORTANT:
open leiden_annotation_sparse.m in MATLAB and paste in line 24
i.e. env_bin = 'F:\Anaconda\envs\leiden_clustering\python.exe';
(Make sure to leave '\python.exe' at the end for windows machine)
Save file
 
 
Usage:
sce = leiden_annotation_sparse(sce,'mouse','knn');
 
BUG: IT WILL FAIL ONCE, RUN AGAIN AND WILL WORK AFTER PYTHON ENVIRONMENT IS SET
 
 
