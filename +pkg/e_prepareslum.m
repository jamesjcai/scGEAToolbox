function e_prepareslum(sce)

if nargin<1, sce=[]; end

answer=questdlg('Select a folder to save the outupt files. Continue?','');
if ~strcmp(answer,'Yes'), return; end

outdir = uigetdir;
if ~isfolder(outdir), return; end

filesaved = fullfile(outdir, 'run_tsne.m');
fid=fopen(filesaved,'w');
fprintf(fid, 'addpath(''../scGEAToolbox'');\n');
fprintf(fid, 'load("clean_data.mat","sce");\n');
fprintf(fid, 'sce=sce.embedcells(''tsne3d'',true,true,2);\n');
fprintf(fid, 'save("clean_data_1.mat","sce",''-v7.3'');\n');
fclose(fid);


filesaved = fullfile(outdir, 'sbatch_run_tsne.slum');
fid=fopen(filesaved,'w');
fprintf(fid,'#!/bin/bash\n');
fprintf(fid,'#SBATCH --export=NONE                #Do not propagate environment\n');
fprintf(fid,'#SBATCH --get-user-env=L             #Replicate login environment\n');
fprintf(fid,'#SBATCH --mail-user=jcai@tamu.edu\n');  
fprintf(fid,'#SBATCH --mail-type=END\n');
fprintf(fid,'#SBATCH --mail-type=FAIL\n');
fprintf(fid,'#SBATCH --job-name=serial_matlab\n');
fprintf(fid,'#SBATCH --time=2:30:00               #Set the wall clock limit to 2hr and 30min\n');
fprintf(fid,'#SBATCH --nodes=1                    #Request 1 node\n');
fprintf(fid,'#SBATCH --ntasks-per-node=8          #Request 8 tasks/cores per node\n');
fprintf(fid,'#SBATCH --mem=32G                    #Request 32GB per node\n');
fprintf(fid,'#SBATCH --output=matlabLogOut.%%j\n');
fprintf(fid,'cd /scratch/user/jcai/matlab/\n');
fprintf(fid,'module purge\n');
fprintf(fid,'module load Matlab/R2022a\n');
fprintf(fid,'matlab -nodisplay -nosplash -nodesktop -batch "run_tsne"\n');
fclose(fid);

if ~isempty(sce)
    filesaved = fullfile(outdir, 'clean_data.mat');
    save(filesaved,"sce");
end