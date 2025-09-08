function [T] = py_scTenifoldXct(sce_ori, celltype1, celltype2, twosided, ...
                                wkdir, isdebug, ...
                                prepare_input_only, parentfig)

T = [];
if nargin < 8, parentfig = []; end
if nargin < 7, prepare_input_only = false; end
if nargin < 6, isdebug = true; end
if nargin < 5, wkdir = []; end
if nargin < 4, twosided = true; end
if nargin < 3, error('Usage: [T] = py_scTenifoldXct(sce, celltype1, celltype2)'); end

sce = copy(sce_ori);

oldpth = pwd();
pw1 = fileparts(mfilename('fullpath'));
codepth = fullfile(pw1, '..', 'external', 'py_scTenifoldXct');

if isempty(wkdir) || ~isfolder(wkdir)
    cd(codepth);
else
    disp('Using working directory provided.');
    cd(wkdir);
end



if ~prepare_input_only
    fw = gui.myWaitbar(parentfig);
    gui.myWaitbar(parentfig, fw, false, [], ...
        'Checking Python environment...');
    
    x = pyenv;
    try
        pkg.i_add_conda_python_path;
        % pyenv('Version', x.Executable, 'ExecutionMode', 'OutOfProcess');
    catch ME
        warning(ME.message);
    end

    codefullpath = fullfile(codepth,'require.py');
    %cmdlinestr = sprintf('"%s" "%s%srequire.py"', ...
    %    x.Executable, codepth, filesep);
    cmdlinestr = sprintf('"%s" "%s"', x.Executable, codefullpath);
    
    disp(cmdlinestr)
    [status, cmdout] = system(cmdlinestr, '-echo');
    if status ~= 0
            
        if isvalid(fw), gui.myWaitbar(parentfig, fw, true); end
        % gui.myErrordlg(parentfig, sprintf('%s', cmdout));
        a = sprintf("%s.", cmdout);
        if strcmp('Yes', gui.myQuestdlg(parentfig, a+" Continue with script.py preparation?"))
            prepare_input_only = true;
        else
            cd(oldpth);
            return;
        end
        % error('Python scTenifoldXct has not been installed properly.');
    end
    if isvalid(fw)
        gui.myWaitbar(parentfig, fw, false, [], 'Checking Python environment is complete');
        pause(1);
        close(fw);
    end
end

tmpfilelist = {'X.mat', 'X.txt', 'g.txt', 'c.txt', 'output.txt', ...
    'output1.txt', 'output2.txt', ...
    'gene_name_Source.tsv', 'gene_name_Target.tsv', ...
    'pcnet_Source.npz', 'pcnet_Target.npz', ...
    'A1.mat', 'A2.mat', 'pcnet_Source.mat', 'pcnet_Target.mat'};

if ~isdebug, pkg.i_deletefiles(tmpfilelist); end

% load(fullfile(pw1,'..','assets','Ligand_Receptor','Ligand_Receptor.mat'), ...
%     'ligand','receptor');
% validg=unique([ligand receptor]);
% [y]=ismember(upper(sce.g),validg);
% X=sce.X(y,:);
% g=sce.g(y);
% writematrix(sce.X,'X.txt');

idx = sce.c_cell_type_tx == celltype1 | sce.c_cell_type_tx == celltype2;
sce = sce.selectcells(idx);  % OK
sce.c_batch_id = sce.c_cell_type_tx;
sce.c_batch_id(sce.c_cell_type_tx == celltype1) = "Source";
sce.c_batch_id(sce.c_cell_type_tx == celltype2) = "Target";
% sce=sce.qcfilter;


if issparse(sce.X)
    X = single(full(sce.X)); 
else
    X = single(sce.X);
end
save('X.mat', '-v7.3', 'X', 'twosided');
writematrix(sce.g, 'g.txt');
writematrix(sce.c_batch_id, 'c.txt');
disp('Input X g c written.');

t = table(sce.g, sce.g, 'VariableNames', {' ', 'gene_name'});
writetable(t, 'gene_name_Source.tsv', 'filetype', 'text', 'Delimiter', '\t');
writetable(t, 'gene_name_Target.tsv', 'filetype', 'text', 'Delimiter', '\t');
disp('Input gene_names written.');

useexist = false;
if exist("pcnet_Source.mat", 'file') && exist("pcnet_Target.mat", 'file')    
    answer = gui.i_questdlgtimer(10, ...
        'pcnet_Source.mat and pcnet_Target.mat files existing. Use them?', ...
         '', 'Yes', ...
         'No', ...
         'Cancel', 'No');
    switch answer
        case 'Yes'
            useexist = true;
        case 'No'
            useexist = false;
        case 'Cancel'
            return;
        otherwise
            return;
    end
end

if ~useexist
    fw = gui.myWaitbar(parentfig);
    gui.myWaitbar(parentfig, fw, false, [], 'Step 1 of 3: Building pcnet_Source network...');
    disp('Building pcnet_Source network...');
    A1 = sc_pcnetpar(sce.X(:, sce.c_cell_type_tx == celltype1));
    A1 = A1 ./ max(abs(A1(:)));
    A = ten.e_filtadjc(A1, 0.75, false);
    save('pcnet_Source.mat', 'A', '-v7.3');
    disp('pcnet_Source.mat saved.');
    if isvalid(fw), gui.myWaitbar(parentfig, fw, false, [], 'Building pcnet_Source is complete'); end
end
if ~useexist
    gui.myWaitbar(parentfig, fw, false, [], 'Step 2 of 3: Building pcnet_Target network...');
    disp('Building pcnet_Target network...')
    A2 = sc_pcnetpar(sce.X(:, sce.c_cell_type_tx == celltype2));
    A2 = A2 ./ max(abs(A2(:)));
    A = ten.e_filtadjc(A2, 0.75, false);
    save('pcnet_Target.mat', 'A', '-v7.3');
    disp('pcnet_Target network saved.')
    if isvalid(fw), gui.myWaitbar(parentfig, fw, false, [], 'Building pcnet_Target is complete'); end
end

if twosided
    twosidedtag = 1;
else
    twosidedtag = 0;
end

if ~prepare_input_only
    gui.myWaitbar(parentfig, fw, false, [], 'Step 3 of 3: Running scTenifoldXct.py...');
else
    gui.myWaitbar(parentfig, fw, false, [], 'Step 3 of 3: Finishing input preparation...');
end

codefullpath = fullfile(codepth,'script.py');
pkg.i_addwd2script(codefullpath, wkdir, 'python');

if ~prepare_input_only
    cmdlinestr = sprintf('"%s" "%s" %d', x.Executable, codefullpath, twosidedtag);
    disp(cmdlinestr)
    [status] = system(cmdlinestr, '-echo');
end
% https://www.mathworks.com/matlabcentral/answers/334076-why-does-externally-called-exe-using-the-system-command-freeze-on-the-third-call
if isvalid(fw)
    if prepare_input_only
        gui.myWaitbar(parentfig, fw, false, [], 'Input preparation is complete.');
    else
        gui.myWaitbar(parentfig, fw, false, [], 'Running scTenifoldXct.py is complete.');
    end
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
if ~prepare_input_only
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
end
%end

if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
cd(oldpth);

if isvalid(fw)
    gui.myWaitbar(parentfig, fw);
end

end
