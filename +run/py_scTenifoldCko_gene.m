function [T] = py_scTenifoldCko_gene(sce_ori, celltype1, celltype2, targetg, ...
                                targettype, wkdir, ...
                                isdebug, prepare_input_only)

T = [];
if nargin < 8, prepare_input_only = false; end
if nargin < 7, isdebug = true; end
if nargin < 6, wkdir = []; end
if nargin < 5 || isempty(targettype), targettype=sprintf('%s+%s', celltype1, celltype2); end
if nargin < 4 || isempty(targetg), targetg = sce_ori.g(1); end
if nargin < 3, error('Usage: [T] = py_scTenifoldCko(sce, celltype1, celltype2, targetg)'); end

twosided = true;

sce = copy(sce_ori);
sce1 = copy(sce);
sce2 = copy(sce);

if targettype==sprintf('%s+%s', celltype1, celltype2)
    sce2.X(sce2.g==targetg, :)=0;
elseif targettype==celltype1
    sce2.X(sce2.g==targetg, sce2.c_cell_type_tx==celltype1)=0;
elseif targettype==celltype2
    sce2.X(sce2.g==targetg, sce2.c_cell_type_tx==celltype2)=0;
else
    error('py_scTenifoldCko error.');
end

%[Tcell, iscomplete] = py_scTenifoldXct2(sce1, sce2, celltype1, celltype2, ...
%                         twosided, wkdir, isdebug, prepare_input_only);


% -----------------

oldpth = pwd();
pw1 = fileparts(mfilename('fullpath'));
codepth = fullfile(pw1, '..', 'external', 'py_scTenifoldCko');

if isempty(wkdir) || ~isfolder(wkdir)
    % cd(codepth);
    wkdir=tempdir;
    cd(wkdir);
else
    disp('Using working directory provided.');
    cd(wkdir);
end

if ~prepare_input_only
 x = pyenv;
    %{
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

    if isvalid(fw)
        gui.gui_waitbar(fw, [], 'Checking Python environment is complete');
    end
    %}
end

    
    tmpfilelist = {'X1.mat', 'X2.mat', 'g1.txt', 'c1.txt', 'g2.txt', 'c2.txt', 'output.txt', ...
        '1/gene_name_Source.tsv', '1/gene_name_Target.tsv', ...
        '2/gene_name_Source.tsv', '2/gene_name_Target.tsv', ...
        '1/pcnet_Source.mat', '1/pcnet_Target.mat', ...
        '2/pcnet_Source.mat', '2/pcnet_Target.mat'};
    
    if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
    
    % load(fullfile(pw1,'..','assets','Ligand_Receptor','Ligand_Receptor.mat'), ...
    %     'ligand','receptor');
    % validg=unique([ligand receptor]);
    % [y]=ismember(upper(sce.g),validg);
    % X=sce.X(y,:);
    % g=sce.g(y);
    % writematrix(sce.X,'X.txt');
    

    in_prepareX(sce1, 1);
    in_prepareX(sce2, 2);

    fw = gui.gui_waitbar([], [], 'Step 2 of 4: Building S1 networks...');
    %try
        in_prepareA12(sce1, targetg);
    % catch ME
    %     if isvalid(fw)
    %         gui.gui_waitbar(fw, [], 'Building S1 networks is incomplete');
    %     end
    %     errordlg(ME.message);
    %     return;
    % end
    gui.gui_waitbar(fw, [], 'Building S1 networks is complete');

    fw = gui.gui_waitbar([], [], 'Step 3 of 4: Building S2 networks...');
    % try
    %     in_prepareA(sce2, 2);
    % catch ME
    %     if isvalid(fw)
    %         gui.gui_waitbar(fw, [], 'Building S2 network is incomplete');
    %     end
    %     errordlg(ME.message);
    %     return;
    % end
    pause(3);
    gui.gui_waitbar(fw, [], 'Building S2 network is complete');

    fw = gui.gui_waitbar([], [], 'Step 4 of 4: Running scTenifoldXct.py...');

codefullpath = fullfile(codepth,'script_gene.py');
pkg.i_addwd2script(codefullpath, wkdir, 'python');


if ~prepare_input_only
    twosidedtag = 2;
    cmdlinestr = sprintf('"%s" "%s" %d', x.Executable, codefullpath, twosidedtag);
    disp(cmdlinestr)
    try
        % [status] = system(cmdlinestr, '-echo');
        [status] = system(cmdlinestr);
        % https://www.mathworks.com/matlabcentral/answers/334076-why-does-externally-called-exe-using-the-system-command-freeze-on-the-third-call
    catch ME
        if isvalid(fw)
            gui.gui_waitbar(fw, [], 'Running scTenifoldCko.py is incomplete.');
        end
        errordlg(ME.message);
        return;
    end
end
    if isvalid(fw)
        if prepare_input_only
            gui.gui_waitbar(fw, [], 'Input preparation is complete.');
        else
            gui.gui_waitbar(fw, [], 'Running scTenifoldCko.py is complete.');
        end
    end    

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
        error('scTenifoldCko runtime error.');
    end
    end

    % if status == 0 && exist('output.txt', 'file')
    %     T = readtable('output.txt');
    %     iscomplete = true;
    % end
    if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
    cd(oldpth);


% --------------------------------------------------
% --------------------------------------------------
% --------------------------------------------------

    function in_prepareX(sce, id)
        if ~exist(sprintf('%d', id), 'dir')
            mkdir(sprintf('%d', id));
        end
        idx = sce.c_cell_type_tx == celltype1 | sce.c_cell_type_tx == celltype2;
        sce = sce.selectcells(idx); % OK
        sce.c_batch_id = sce.c_cell_type_tx;
        sce.c_batch_id(sce.c_cell_type_tx == celltype1) = "Source";
        sce.c_batch_id(sce.c_cell_type_tx == celltype2) = "Target";
        % sce=sce.qcfilter;
        if issparse(sce.X)
            X = single(full(sce.X)); 
        else
            X = single(sce.X);
        end
        twosided = true;
        save(sprintf('X%d.mat', id), '-v7.3', 'X', 'twosided');
        writematrix(sce.g, sprintf('g%d.txt', id));
        writematrix(sce.c_batch_id, sprintf('c%d.txt', id));
        fprintf('Input X%d g%d c%d written.\n', id, id, id);
        t = table(sce.g, sce.g, 'VariableNames', {' ', 'gene_name'});
        writetable(t, sprintf('%d/gene_name_Source.tsv', id), ...
            'filetype', 'text', 'Delimiter', '\t');
        writetable(t, sprintf('%d/gene_name_Target.tsv', id), ...
            'filetype', 'text', 'Delimiter', '\t');
        disp('Input gene_names written.');
    end

    function in_prepareA12(sce, targetg)

        idx = find(sce.g==targetg);
        assert(isscalar(idx));

        disp('Building A1 network...')
        A1 = sc_pcnetpar(sce.X(:, sce.c_cell_type_tx == celltype1));
        disp('A1 network built.')
        A1 = A1 ./ max(abs(A1(:)));
        % A=0.5*(A1+A1.');
        A = ten.e_filtadjc(A1, 0.75, false);
        save(sprintf('%d/pcnet_Source.mat', 1), 'A', '-v7.3');

        if contains(targettype, celltype1)
            fprintf('\nKO gene in %s.\n', celltype1);
            if nnz(A(idx, :) ~= 0) < 10
                warning('KO gene (%s) has no link or too few links with other genes.', ...
                    targetg);
            end
            A(idx,:)=0;
        end
        save(sprintf('%d/pcnet_Source.mat', 2), 'A', '-v7.3');


        disp('Building A2 network...');
        A2 = sc_pcnetpar(sce.X(:, sce.c_cell_type_tx == celltype2));
        disp('A2 network built.');
        A2 = A2 ./ max(abs(A2(:)));
        A = ten.e_filtadjc(A2, 0.75, false);
        save(sprintf('%d/pcnet_Target.mat', 1), 'A', '-v7.3');

        if contains(targettype, celltype2)
            fprintf('\nKO gene in %s.\n', celltype2);
            if nnz(A(idx, :) ~= 0) < 10
                warning('KO gene (%s) has no link or too few links with other genes.', ...
                    targetg);
            end
            A(idx,:)=0;
        end
        save(sprintf('%d/pcnet_Target.mat', 2), 'A', '-v7.3');
    end
     
    % function in_prepareA(sce, id)
    %     disp('Building A1 network...')
    %     A1 = sc_pcnetpar(sce.X(:, sce.c_cell_type_tx == celltype1));
    %     disp('A1 network built.')
    %     A1 = A1 ./ max(abs(A1(:)));
    %     % A=0.5*(A1+A1.');
    %     A = ten.e_filtadjc(A1, 0.75, false);
    %     save(sprintf('%d/pcnet_Source.mat', id), 'A', '-v7.3');
    % 
    %     disp('Building A2 network...');
    %     A2 = sc_pcnetpar(sce.X(:, sce.c_cell_type_tx == celltype2));
    %     disp('A2 network built.');
    %     A2 = A2 ./ max(abs(A2(:)));
    %     % A=0.5*(A2+A2.');
    %     A = ten.e_filtadjc(A2, 0.75, false);
    %     save(sprintf('%d/pcnet_Target.mat', id), 'A', '-v7.3');
    % end

end