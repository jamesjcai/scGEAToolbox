function [done] = i_setqtmservice(src, ~)

[parentfig, ~] = gui.gui_getfigsce(src);
done = false;

installedPackages = matlabshared.supportpkg.getInstalled;
isQuantumInstalled = any(strcmp({installedPackages.Name}, 'MATLAB Support Package for Quantum Computing'));

if isQuantumInstalled
    % disp('Quantum Computing Support Package is already installed.');
    gui.myHelpdlg(parentfig, "Quantum Computing Support Package is already installed.");

else

    %url = 'https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/53c14e1e-1831-4a74-b0ee-e306c3cbf83d/7e5d5d80-5c7c-42ef-b5d9-3944e9549d0b/packages/mlpkginstall';
    %tempZip = fullfile(tempdir, 'MLSupportPackageForQuantumComputing.mlpkginstall');
    %websave(tempZip, url);
    %matlab.addons.install(tempZip);
    if strcmpi('Yes', gui.myQuestdlg(parentfig, 'Quantum Computing Support Package is not installed. Visit the MATLAB Support Package for Quantum Computing page to download the installer?'))
        web('https://www.mathworks.com/matlabcentral/fileexchange/125425-matlab-support-package-for-quantum-computing');
    end
end


%{
preftagname = 'qtmserviceapi';
if ~ispref('scgeatoolbox', preftagname)
    % answer = gui.myQuestdlg(parentfig, 'LLM model has not been set up. Set it up?');
    % if ~strcmp(answer, 'Yes'), return; end
    % [done] = ix_setwdpath(pathdefult);
else
    s = getpref('scgeatoolbox', preftagname);
    answer1 = gui.myQuestdlg(parentfig, sprintf('%s', s), ...
        'Quantum Service API', ...
        {'Use this', 'Use another', 'Cancel'}, 'Use this');
    if isempty(answer1), return; end
    switch answer1
        case 'Use this'
            done = true;
            return;
        case 'Use another'
            % [done] = ix_setwdpath(s);
        otherwise
            return;
    end
end

listItems = {'IBM Quantum', 'Amazon Braket', 'Google Quantum',...
    'Microsoft Azure Quantum', 'Xanadu PennyLane'};

if gui.i_isuifig(parentfig)
    [selectedIndex, ok] = gui.myListdlg(parentfig, listItems, ...
            'Select a Quantum Service API:', listItems(1));
else
    [selectedIndex, ok] = listdlg('PromptString', 'Select Quantum Service API:', ...
                          'SelectionMode', 'single', ...
                          'ListString', listItems, ...
                          'ListSize', [220 300], ...
                          'InitialValue', 1);
end

if ok
    selectedProvider = listItems{selectedIndex};
    % fprintf('You selected: %s\n', selectedProvider);
    switch selectedProvider
        case 'Ollama'
            done = true;
        otherwise
            gui.myWarndlg(parentfig, ...
                sprintf(['The function supporting %s API is ' ...
                'under development.'], ...
                selectedProvider));
            return;
    end
else
    return;
end

if done
     gui.myHelpdlg(parentfig, "Quantum Service API has been set successfully.");
end
%}