function sc_llm_enrichr2word(selpath, parentfig)

if nargin<2, parentfig = []; end
if nargin < 1
    selpath = uigetdir;
end
if ~isempty(parentfig) && pkg.i_isvalid(parentfig) && parentfig.Visible == "on"
    figure(parentfig);
    cleanupObj = onCleanup(@() gui.i_raisefig(parentfig));
end
if isempty(selpath) || isequal(selpath, 0), return; end
if ~isfolder(selpath), return; end



files = dir(fullfile(selpath, '*DE_*.xlsx'));
fileNames1 = string({files(~[files.isdir]).name});
files = dir(fullfile(selpath, '*DV_*.xlsx'));
fileNames2 = string({files(~[files.isdir]).name});

listItems = [fileNames1'; fileNames2'];

if isempty(listItems)
    gui.myHelpdlg(parentfig, ...
        'No DE/DV Excel files with Enrichr results were found in the selected folder.');
    return;
end

if gui.i_isuifig(parentfig)
    [selectedIndex, ok] = gui.myListdlg(parentfig, listItems, ...
            'Select Excel Files:', listItems, true);
else
    [selectedIndex, ok] = listdlg('PromptString', 'Select Excel Files:', ...
                          'SelectionMode', 'multiple', ...
                          'ListString', listItems, ...
                          'ListSize', [260 300], ...
                          'InitialValue', 1:numel(listItems));
end

if ok
    selectedfiles = listItems(selectedIndex);
else
    return;
end

fw = gui.myWaitbar(parentfig);
closeFw = onCleanup(@() gui.myWaitbar(parentfig, fw, true));
nwritten = 0;
failed = strings(0, 1);

for k = 1:length(selectedfiles)
    gui.myWaitbar(parentfig, fw, false, '', ...
        sprintf('%s', selectedfiles(k)), ...
        (k-0.5)/length(selectedfiles));
    infile = fullfile(selpath, selectedfiles(k));
    % [TbpUpEnrichr, TmfUpEnrichr, ...
    %     TbpDnEnrichr, TmfDnEnrichr] = in_gettables(infile);

    [~, wordfilename] = fileparts(selectedfiles(k));
    % Per file: an unreadable workbook or a failed LLM call used to stop
    % the loop with the bar open, and a silent failure told no one.
    try
        [TbpUpEnrichr, TmfUpEnrichr, ...
            TbpDnEnrichr, TmfDnEnrichr] = pkg.in_XLSX2DETable(infile);
        done = llm.e_DETableSummary(TbpUpEnrichr, ...
            TmfUpEnrichr, TbpDnEnrichr, ...
            TmfDnEnrichr, wordfilename, selpath);
    catch ME
        fprintf('Report for %s failed: %s\n', selectedfiles(k), ME.message);
        done = false;
    end
    if done
        nwritten = nwritten + 1;
    else
        failed(end+1) = selectedfiles(k); %#ok<AGROW>
    end
end
gui.myWaitbar(parentfig, fw);
gui.i_reportllmword(parentfig, nwritten, failed, selpath);

end
