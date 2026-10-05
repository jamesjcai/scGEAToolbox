function [insce, names] = i_loadscefiles(parentfig, pathname, fnames)
% I_LOADSCEFILES - Read the SCE out of each of several .mat files.
%
%   [insce, names] = gui.i_loadscefiles(parentfig, pathname, fnames)
%
% FNAMES is a file name or a cell array of them, as UIGETFILE returns, all in
% folder PATHNAME. INSCE is a cell array with one SingleCellExperiment per
% file and NAMES the file names without extension, as a string array.
%
% A file is used when it holds a variable named sce, or exactly one SCE under
% any other name. Anything else is an error that names the file, rather than
% what a bare LOAD(file, 'sce') does: nothing at all for the first file, and
% for each later one, quietly reuse the SCE read from the file before it.
%
% Each SCE records the file it came from in its metadata. Returns empty
% outputs, having shown the error, when a file cannot be used.
%
% see also: gui.i_mergesces, gui.sc_openscedlg

fnames = cellstr(fnames);
insce = cell(1, numel(fnames));
names = strings(1, numel(fnames));

fw = gui.myWaitbar(parentfig);
try
    for k = 1:numel(fnames)
        scefile = fullfile(pathname, fnames{k});
        % LOAD reports no progress, so the bar can only count files already
        % read. With one file that count says nothing, and no fraction
        % leaves the bar spinning.
        if isscalar(fnames)
            frac = [];
        else
            frac = (k - 1)/numel(fnames);
        end
        gui.myWaitbar(parentfig, fw, false, '', ...
            sprintf('Loading %s...', fnames{k}), frac);
        info = whos('-file', scefile);
        scevars = {info(strcmp({info.class}, 'SingleCellExperiment')).name};
        if ismember('sce', scevars)
            varname = 'sce';
        elseif isscalar(scevars)
            varname = scevars{1};
        elseif isempty(scevars)
            error('i_loadscefiles:NoSCEInFile', ...
                '%s does not contain an SCE (SingleCellExperiment) variable.', ...
                fnames{k});
        else
            error('i_loadscefiles:SeveralSCEsInFile', ...
                ['%s holds %d SCE variables, none named sce. Save the ' ...
                'dataset to use as a variable named sce.'], ...
                fnames{k}, numel(scevars));
        end
        data = load(scefile, varname);
        insce{k} = data.(varname).appendmetainfo(sprintf("Source: %s", scefile));
        [~, names(k)] = fileparts(fnames{k});
    end
catch ME
    gui.myWaitbar(parentfig, fw, true);
    gui.myErrordlg(parentfig, ME.message, ME.identifier);
    insce = {};
    names = strings(0);
    return;
end
gui.myWaitbar(parentfig, fw);
end
