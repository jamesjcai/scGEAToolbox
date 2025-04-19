function [Tbp, Tmf] = py_GSEApy_enr(genelist, backgroundlist, wkdir, ...
    showbarplot, showprogress, isdebug)

if nargin < 6, isdebug = false; end
if nargin < 5, showprogress = false; end
if nargin < 4, showbarplot = false; end
if nargin < 3, wkdir = tempdir; end
if nargin < 2, backgroundlist = []; end
if nargin < 1, return; end
Tmf=[]; Tbp = [];

extprogname = 'py_GSEApy_enr';
if isempty(wkdir)
    preftagname = 'externalwrkpath';
    [wkdir] = gui.gui_setprgmwkdir(extprogname, preftagname);
    if isempty(wkdir), return; end
end

    oldpth = pwd();
    pw1 = fileparts(mfilename('fullpath'));
    codepth = fullfile(pw1, 'external', extprogname);
    if isempty(wkdir) || ~isfolder(wkdir)
        cd(codepth);
    else
        % disp('Using working directory provided.');
        cd(wkdir);
    end

    if showprogress
        fw = gui.gui_waitbar([], [], 'Checking Python environment...');
    end

    x = pyenv;
    try
        pkg.i_add_conda_python_path;
    catch
    
    end
    
    codefullpath = fullfile(codepth,'require.py');
    cmdlinestr = sprintf('"%s" "%s"', x.Executable, codefullpath);
    disp(cmdlinestr)
    [status, cmdout] = system(cmdlinestr, '-echo');
    if status ~= 0
        cd(oldpth);
        if showprogress && isvalid(fw)
            gui.gui_waitbar(fw, true);
        end
        error(cmdout);
    end


tmpfilelist = {'input.txt', 'background.txt',...
    'GO_Biological_Process_2025.Human.enrichr.reports.txt',...
    'GO_Molecular_Function_2025.Human.enrichr.reports.txt',...
    'GO_Biological_Process_2025.Human.enrichr.reports.png',...
    'GO_Molecular_Function_2025.Human.enrichr.reports.png'};
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end

usebktag = 1;
writematrix(genelist, 'input.txt');
if ~isempty(backgroundlist)
    writematrix(backgroundlist, 'background.txt');
    usebktag = 2;
    disp('Background file written.');
end

if showprogress && isvalid(fw)
    gui.gui_waitbar(fw, [], [], 'Checking Python environment is complete');
    pause(0.5);
    gui.gui_waitbar(fw, [], [], 'Running GSEApy Enrichr...');
end
% fw = gui.gui_waitbar([],[],'Running DataMapPlot...');

    codefullpath = fullfile(codepth,'script.py');
    pkg.i_addwd2script(codefullpath, wkdir, 'python');
    
    cmdlinestr = sprintf('"%s" "%s" %d', x.Executable, codefullpath, usebktag);

    disp(cmdlinestr)
    [status] = system(cmdlinestr, '-echo');

if status == 0 && showprogress && isvalid(fw)
    gui.gui_waitbar(fw, [], 'GSEApy Enrichr is complete');
end
if exist('GO_Biological_Process_2025.Human.enrichr.reports.txt','file')
    Tbp = readtable('GO_Biological_Process_2025.Human.enrichr.reports.txt', ...
            'VariableNamingRule', 'modify');
    [Tbp] = in_trimT(Tbp);
end
if exist('GO_Molecular_Function_2025.Human.enrichr.reports.txt','file')
    Tmf = readtable('GO_Molecular_Function_2025.Human.enrichr.reports.txt', ...
        'VariableNamingRule', 'modify');
    [Tmf] = in_trimT(Tmf);
end

if showbarplot
    if exist('GO_Biological_Process_2025.Human.enrichr.reports.png','file')
        figure; imshow(imread('GO_Biological_Process_2025.Human.enrichr.reports.png'))
    end
    if exist('GO_Molecular_Function_2025.Human.enrichr.reports.png','file')
        figure; imshow(imread('GO_Molecular_Function_2025.Human.enrichr.reports.png'))
    end
end

if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
cd(oldpth);
end


function [T] = in_trimT(T)
    %try
        if ~isempty(T)
            if ismember('Overlap', T.Properties.VariableNames) && ismember('P_value', T.Properties.VariableNames)
                a = string(T.Overlap);
                b = extractBefore(a, cell2mat(strfind(a,'/')));
                idx = double(b)>=5 & T.P_value < 0.05;
                T = T(idx,:);
            elseif ismember('Genes', T.Properties.VariableNames) && ismember('P_value', T.Properties.VariableNames)
                idx = cellfun(@numel,strfind(T.Genes,';'))+1 >= 5 & T.P_value < 0.05;
                T = T(idx,:);
            end
        end
    %catch
    %end
end