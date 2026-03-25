function i_installquantumpackage(src, ~)

[parentfig, ~] = gui.gui_getfigsce(src);

if ~(ismcc || isdeployed)
    %#exclude matlabshared.supportpkg.getInstalled
    installedPackages = matlabshared.supportpkg.getInstalled;
else
    installedPackages = [];
end

if isempty(installedPackages)
    isQuantumInstalled = false;
else
    isQuantumInstalled = any(strcmp({installedPackages.Name}, 'MATLAB Support Package for Quantum Computing'));
end

if isQuantumInstalled
    gui.myHelpdlg(parentfig, "Quantum Computing Support Package is already installed.");
else
    if strcmpi('Yes', gui.myQuestdlg(parentfig, 'Quantum Computing Support Package is not installed. Visit the MATLAB Support Package for Quantum Computing page to download the installer?'))
        web('https://www.mathworks.com/matlabcentral/fileexchange/125425-matlab-support-package-for-quantum-computing');
    end
end
