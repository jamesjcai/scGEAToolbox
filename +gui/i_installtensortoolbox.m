function i_installtensortoolbox(src, ~)

[parentfig] = gui.gui_getfigsce(src);

if exist(['@tensor', filesep, 'tensor.m'], 'file') == 2
    gui.myHelpdlg(parentfig, 'Tensor Toolbox is already installed.');
    return;
end

answer = gui.myQuestdlg(parentfig, ...
    ['Tensor Toolbox for MATLAB (Sandia Labs) is required but not installed. ' ...
     'Download and install automatically (requires internet)?'], ...
    'Tensor Toolbox', ...
    {'Download', 'Visit Website', 'Cancel'}, 'Download');

switch answer
    case 'Download'
        installBase = fullfile(prefdir, 'scgeatoolbox_addons');
        if ~isfolder(installBase), mkdir(installBase); end
        url = 'https://github.com/sandialabs/tensor_toolbox/archive/refs/heads/master.zip';
        tmpZip = fullfile(tempdir, 'tensor_toolbox.zip');
        fw = gui.myWaitbar(parentfig);
        try
            websave(tmpZip, url);
            unzip(tmpZip, installBase);
            gui.myWaitbar(parentfig, fw);
            d = dir(fullfile(installBase, 'tensor_toolbox*'));
            d = d([d.isdir]);
            if isempty(d)
                error('Extraction failed: tensor_toolbox folder not found.');
            end
            pth = fullfile(installBase, d(1).name);
            addpath(pth);
            setpref('scgeatoolbox', 'tensor_toolbox_path', pth);
            gui.myHelpdlg(parentfig, 'Tensor Toolbox installed successfully.');
        catch ME
            gui.myWaitbar(parentfig, fw);
            gui.myErrordlg(parentfig, ['Download failed: ' ME.message]);
        end
    case 'Visit Website'
        web('https://www.tensortoolbox.org', '-browser');
end
end
