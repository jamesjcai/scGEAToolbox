clc;
disp('% 1/5 - Change current directory >>cdgea;')
fprintf('Press any key to continue...')
pause
cdgea;
fprintf('...DONE.\n');
pause(2)
clc

disp('% 2/5 - Load example data 1')
disp('>>load example_data/testXgs.mat X g s')
disp('% X=expression matrix; g=gene list; s=embedding') 
fprintf('Press any key to continue...')
pause
load example_data/testXgs.mat X g s
fprintf('...DONE.\n\n');
pause(2)
clc

disp('% 3/5 - Load example data 2')
disp('>>load example_data/testSce.mat sce')
disp('% sce=class SingleCellExperiment') 
fprintf('Press any key to continue...')
pause
load example_data/testSce sce
fprintf('...DONE.\n');
pause(2)
clc

disp('% 4/5 - Run SC_SCATTER program >>sc_scatter(X,g,s);')
disp('% 5/5 - Run SC_SCATTER_SCE program >>sc_scatter_sce(sce);')

fprintf('Press any key to continue...')
pause
sc_scatter(X,g,s);
%fprintf('...DONE.\n\n');


%fprintf('Press any key to continue...')
pause(3)
f=sc_scatter_sce(sce);
f.Position(1)=f.Position(1)+50;
f.Position(2)=f.Position(2)-50;
fprintf('...DONE.\n\n');

pause(3)
gui.sc_multiembeddings