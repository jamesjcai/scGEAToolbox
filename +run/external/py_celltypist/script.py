import os
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
# https://www.blogforbrains.com/blog/2014/9/6/loading-matlab-mat-data-in-python
# import scipy.io as spio

import celltypist
input_file = 'input.csv'
predictions = celltypist.annotate(input_file, majority_voting = True, model='Developing_Mouse_Brain.pkl')
predictions.predicted_labels.to_csv('output.csv')

# https://www.science.org/doi/epdf/10.1126/science.abl5197
#  Downloading model [1/12]: Immune_All_Low.pkl
#  Downloading model [2/12]: Immune_All_High.pkl
#  Downloading model [3/12]: Adult_Mouse_Gut.pkl
#  Downloading model [4/12]: COVID19_Immune_Landscape.pkl
#  Downloading model [5/12]: Cells_Fetal_Lung.pkl
#  Downloading model [6/12]: Cells_Intestinal_Tract.pkl
#  Downloading model [7/12]: Cells_Lung_Airway.pkl
#  Downloading model [8/12]: Developing_Mouse_Brain.pkl
#  Downloading model [9/12]: Healthy_COVID19_PBMC.pkl
#  Downloading model [10/12]: Human_Lung_Atlas.pkl
#  Downloading model [11/12]: Nuclei_Lung_Airway.pkl
#  Downloading model [12/12]: Pan_Fetal_Human.pkl
