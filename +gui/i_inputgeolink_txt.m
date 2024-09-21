function [X, genelist, celllist, ftdone, answer1] = i_inputgeolink_txt

X = [];
genelist = [];
celllist = [];
ftdone = false;
answer1 = [];

prompt = {'Enter link to counts.txt.gz or counts.csv.gz:'};
dlgtitle = 'Input Download Links';
dims = [1, 100];
definput = {'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5350nnn/GSM5350808/suppl/GSM5350808_Fibroblast_young_1wk_Saline_counts.csv.gz'};
answer = inputdlg(prompt, dlgtitle, dims, definput);
if isempty(answer), return; end

if ~(ismcc || isdeployed)
    %#exclude urldecode
    answer1 = urldecode(answer{1});
else
    answer1 = pkg.urldecoding(answer{1});
end

if ~isempty(answer1)
    fw = gui.gui_waitbar;
    try
        tmpd = tempdir;
        if strcmpi(answer1(end-2:end), '.gz')
            fprintf('gunzip(''%s'',''%s'');\n', answer1, tmpd);
            files = gunzip(answer1, tmpd);
        elseif strcmpi(answer1(end-3:end), '.zip')
            fprintf('unzip(''%s'',''%s'');\n', answer1, tmpd);
            files = unzip(answer1, tmpd);
        elseif strcmpi(answer1(end-3:end), '.csv')
            files = websave(tempname, answer1);
        else
            error('File format is not supported.');
        end
        if iscell(files)
            f = files{1};
        else
            f = files;
        end
        f = strrep(f, '?', '_');

        if isempty(f) || ~exist(f, 'file')
            error('File format is not supported.');
        end

        fprintf('[X,g]=sc_readtsvfile(''%s'');\n', f);
        [X, genelist, celllist] = sc_readtsvfile(f);
        ftdone = true;
    catch ME
        gui.gui_waitbar(fw, true);
        errordlg(ME.message);
        return;
    end
    gui.gui_waitbar(fw);
end
end
