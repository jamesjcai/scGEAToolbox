function i_installtensortoolbox(src, ~)

[parentfig] = gui.gui_getfigsce(src);

% TEN.CHECK_TENSOR_TOOLBOX also finds a copy this function saved earlier
% (the tensor_toolbox_path preference) and puts it on the path. Looking at
% the path alone missed it in every new session, and offered to download
% it again.
try
    ten.check_tensor_toolbox;
    gui.myHelpdlg(parentfig, 'Tensor Toolbox is already installed.');
    return;
catch
    % Not installed: offer it below.
end

answer = gui.myQuestdlg(parentfig, ...
    ['Tensor Toolbox for MATLAB (Sandia Labs) is not installed. ' ...
     'Download and install it automatically (requires internet)?'], ...
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
