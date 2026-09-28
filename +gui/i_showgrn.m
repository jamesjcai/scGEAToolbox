function [hFig, fname] = i_showgrn(A, g, method, tag, parentfig)
% I_SHOWGRN  Save a built GRN to the working folder and draw it.
%
%   [hFig, fname] = gui.i_showgrn(A, g, method, tag, parentfig) saves A
%   and g (and METHOD, when not empty) to grn_<TAG>_<timestamp>.mat in the
%   working folder - Setup > Set Working Folder, or the current folder when
%   that is unset - without asking. It then draws the whole network with
%   sc_grnview when it has at most 100 genes, and a module sketch with
%   sc_grnsketch otherwise. The figure's toolbar gains two buttons: Export
%   Network to Workspace, and Show Network File in Folder.
%
%   FNAME is empty when the file could not be written; the network is
%   still drawn and exportable.
%
% See also gui.callback_BuildGeneNetwork, gui.callback_BuildGRNAllGenes,
%   sc_grnview, sc_grnsketch.

maxWholeGenes = 100;
g = string(g(:));

fname = i_savefile(A, g, method, tag, parentfig);

if numel(g) <= maxWholeGenes
    hFig = sc_grnview(A, g, '', parentfig);
else
    hFig = sc_grnsketch(A, g, ParentFig=parentfig);
end
if isempty(hFig) || ~pkg.i_isvalid(hFig), return; end

tb = findall(hFig, 'Type', 'uitoolbar');
if isempty(tb)
    tb = uitoolbar(hFig);
end
tb = tb(1);
if ~(ismcc || isdeployed)
    gui.i_addbutton2fig(tb, 'on', @in_export, 'export.gif', ...
        'Export Network to Workspace');
end
pt = gui.i_addbutton2fig(tb, 'off', @in_locate, i_foldericon(), ...
    'Show Network File in Folder');
if isempty(fname)
    pt.Enable = 'off';
else
    pt.TooltipString = sprintf('Show Network File in Folder (%s)', fname);
end

    function in_export(~, ~)
        export2wsdlg({'Save network to variable named:', ...
            'Save gene list to variable named:'}, {'A', 'g'}, {A, g});
    end

    function in_locate(~, ~)
        if ~isfile(fname)
            gui.myErrordlg(hFig, sprintf(['%s no longer exists. It may ', ...
                'have been moved or deleted.'], fname));
            return;
        end
        i_revealfile(fname);
    end
end


function fname = i_savefile(A, g, method, tag, parentfig)
% Write without asking; a failure is reported in the Command Window only,
% so the network is still drawn. A genome-wide network takes seconds to
% compress into a v7.3 file, so large ones get a waitbar; small ones would
% only flash it (myWaitbar pauses half a second on opening).
minWaitbarGenes = 1000;
folder = getpref('scgeatoolbox', 'externalwrkpath', '');
if isempty(folder) || ~isfolder(folder)
    folder = pwd;
end
fname = fullfile(folder, sprintf('grn_%s_%s.mat', tag, ...
    string(datetime("now", Format="yyyyMMdd_HHmmss"))));
fw = [];
if numel(g) >= minWaitbarGenes
    fw = gui.myWaitbar(parentfig, [], false, 'Saving network...');
end
try
    if isempty(method)
        save(fname, 'A', 'g', '-v7.3');
    else
        save(fname, 'A', 'g', 'method', '-v7.3');
    end
    if ~isempty(fw), gui.myWaitbar(parentfig, fw); end
    fprintf('The network has been saved in %s\n', fname);
catch ME
    if ~isempty(fw), gui.myWaitbar(parentfig, fw, true); end
    warning('gui:i_showgrn:saveFailed', ...
        'The network could not be saved in %s: %s', folder, ME.message);
    fname = '';
end
end


function i_revealfile(fname)
% Open the containing folder with the file selected where the OS allows it
if ispc
    system(sprintf('explorer.exe /select,"%s"', fname));
elseif ismac
    system(sprintf('open -R "%s"', fname));
else
    system(sprintf('xdg-open "%s" &', fileparts(fname)));
end
end


function img = i_foldericon()
% 16x16 folder glyph; no folder icon ships in assets/Images
bg = 0.94;
img = bg*ones(16, 16, 3);
body = [0.93, 0.74, 0.25];
edge = [0.62, 0.45, 0.08];
tab = false(16);
tab(3:4, 2:7) = true;
box = false(16);
box(5:14, 2:15) = true;
for c = 1:3
    layer = img(:, :, c);
    layer(tab | box) = body(c);
    outline = (box & ~i_erode(box)) | (tab & ~i_erode(tab));
    layer(outline) = edge(c);
    img(:, :, c) = layer;
end
end


function m = i_erode(m)
% 4-neighbour erosion without Image Processing Toolbox
p = false(size(m) + 2);
p(2:end - 1, 2:end - 1) = m;
m = m & p(1:end - 2, 2:end - 1) & p(3:end, 2:end - 1) & ...
    p(2:end - 1, 1:end - 2) & p(2:end - 1, 3:end);
end
