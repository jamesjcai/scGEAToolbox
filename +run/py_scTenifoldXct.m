function [T] = py_scTenifoldXct(sce, celltype1, celltype2, twosided, ...
                                wkdir, isdebug)

if nargin<6, isdebug = true; end

T = [];
if nargin < 5, wkdir = []; end
if nargin < 4, twosided = true; end

oldpth = pwd();
pw1 = fileparts(mfilename('fullpath'));
codepth = fullfile(pw1, 'external', 'py_scTenifoldXct');

if isempty(wkdir) || ~isfolder(wkdir)
    cd(codepth);
else
    disp('Using working directory provided.');
    cd(wkdir);
end



fw = gui.gui_waitbar([], [], 'Checking Python environment...');

x = pyenv;
try
    pkg.i_add_conda_python_path;
catch

end

codefullpath = fullfile(codepth,'require.py');
%cmdlinestr = sprintf('"%s" "%s%srequire.py"', ...
%    x.Executable, codepth, filesep);
cmdlinestr = sprintf('"%s" "%s"', x.Executable, codefullpath);

disp(cmdlinestr)
[status, cmdout] = system(cmdlinestr, '-echo');
if status ~= 0
    cd(oldpth);

    if isvalid(fw)
        gui.gui_waitbar(fw, true);
    end
    %waitfor(errordlg(sprintf('%s',cmdout)));
    error(cmdout);
    %error('Python scTenifoldXct has not been installed properly.');
end


tmpfilelist = {'X.mat', 'X.txt', 'g.txt', 'c.txt', 'output.txt', ...
    'output1.txt', 'output2.txt', ...
    'gene_name_Source.tsv', 'gene_name_Target.tsv', ...
    'pcnet_Source.npz', 'pcnet_Target.npz', ...
    'A1.mat', 'A2.mat', 'pcnet_Source.mat', 'pcnet_Target.mat'};

if ~isdebug, pkg.i_deletefiles(tmpfilelist); end

% load(fullfile(pw1,'..','resources','Ligand_Receptor','Ligand_Receptor.mat'), ...
%     'ligand','receptor');
% validg=unique([ligand receptor]);
% [y]=ismember(upper(sce.g),validg);
% X=sce.X(y,:);
% g=sce.g(y);
% writematrix(sce.X,'X.txt');

idx = sce.c_cell_type_tx == celltype1 | sce.c_cell_type_tx == celltype2;
sce = sce.selectcells(idx);
sce.c_batch_id = sce.c_cell_type_tx;
sce.c_batch_id(sce.c_cell_type_tx == celltype1) = "Source";
sce.c_batch_id(sce.c_cell_type_tx == celltype2) = "Target";
% sce=sce.qcfilter;


if issparse(sce.X)
    X = single(full(sce.X)); 
else
    X = single(sce.X);
end
save('X.mat', '-v7.3', 'X');
writematrix(sce.g, 'g.txt');
writematrix(sce.c_batch_id, 'c.txt');
disp('Input X g c written.');

t = table(sce.g, sce.g, 'VariableNames', {' ', 'gene_name'});
writetable(t, 'gene_name_Source.tsv', 'filetype', 'text', 'Delimiter', '\t');
writetable(t, 'gene_name_Target.tsv', 'filetype', 'text', 'Delimiter', '\t');
disp('Input gene_names written.');

if isvalid(fw)
    gui.gui_waitbar(fw, [], 'Checking Python environment is complete');
end

useexist = false;
if exist("pcnet_Source.mat", 'file')
    answer = gui.i_questdlgtimer(10, ...
        'pcnet_Source.mat existing. Use it?','', 'Yes, use pcnet_Source', 'No, reconstruct pcnet_Source', ...
        'Cancel', 'Yes, use pcnet_Source');
    switch answer
        case 'Yes, use pcnet_Source'
            useexist = true;
        case 'No, reconstruct pcnet_Source'
            useexist = false;
        case 'Cancel'
            return;
        otherwise
            return;
    end
end
if ~useexist
    fw = gui.gui_waitbar([], [], 'Step 1 of 3: Building pcnet\_Source network...');
    disp('Building pcnet_Source network...');
    A1 = sc_pcnetpar(sce.X(:, sce.c_cell_type_tx == celltype1));
    A1 = A1 ./ max(abs(A1(:)));
    A = ten.e_filtadjc(A1, 0.75, false);
    save('pcnet_Source.mat', 'A', '-v7.3');
    disp('pcnet_Source.mat saved.');
end

if isvalid(fw)
    gui.gui_waitbar(fw, [], 'Building pcnet_Source is complete');
end


useexist = false;
if exist("pcnet_Target.mat", 'file')
    answer = gui.i_questdlgtimer(10, ...
        'pcnet\_Target.mat existing. Use it?','', 'Yes, use pcnet_Target', 'No, reconstruct pcnet_Target', ...
        'Cancel', 'Yes, use pcnet_Target');
    switch answer
        case 'Yes, use pcnet_Target'
            useexist = true;
        case 'No, reconstruct pcnet_Target'
            useexist = false;
        case 'Cancel'
            return;
        otherwise
            return;
    end
end
if ~useexist
    fw = gui.gui_waitbar([], [], 'Step 2 of 3: Building pcnet\_Target network...');
    disp('Building pcnet_Target network...')
    A2 = sc_pcnetpar(sce.X(:, sce.c_cell_type_tx == celltype2));
    A2 = A2 ./ max(abs(A2(:)));
    A = ten.e_filtadjc(A2, 0.75, false);
    save('pcnet_Target.mat', 'A', '-v7.3');
    disp('pcnet_Target network saved.')
end

if isvalid(fw)
    gui.gui_waitbar(fw, [], 'Building pcnet\_Target is complete');
end

if twosided
    twosidedtag = 2;
else
    twosidedtag = 1;
end

fw = gui.gui_waitbar([], [], 'Step 3 of 3: Running scTenifoldXct.py...');

codefullpath = fullfile(codepth,'script.py');
pkg.i_addwd2script(codefullpath, wkdir, 'python');
cmdlinestr = sprintf('"%s" "%s" %d', x.Executable, codefullpath, twosidedtag);

% 
% cmdlinestr = sprintf('"%s" "%s%sscript.py" %d', ...
%     x.Executable, codepth, filesep, tag);


disp(cmdlinestr)
[status] = system(cmdlinestr, '-echo');
% https://www.mathworks.com/matlabcentral/answers/334076-why-does-externally-called-exe-using-the-system-command-freeze-on-the-third-call
if isvalid(fw)
    gui.gui_waitbar(fw, [], 'Running scTenifoldXct.py is complete');
end

% rt=java.lang.Runtime.getRuntime();
% pr = rt.exec(cmdlinestr);
% [status]=pr.waitFor();

% if twosided
%     if status==0 && exist('output1.txt','file') && exist('output2.txt','file')
%         T1=readtable('output1.txt');
%         T2=readtable('output2.txt');
%         T={T1,T2};
%     end
% else
if status == 0 && exist('output1.txt', 'file')
    T = readtable('output1.txt');
    if twosided && exist('output2.txt', 'file')
        T2 = readtable('output2.txt');
        T = {T, T2};
    end
else
    if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
    cd(oldpth);
    error('scTenifoldXct runtime error.');
end
%end

if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
cd(oldpth);
end
