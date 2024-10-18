function [outgenelist, outbackgroundlist, enrichrtype] = gui_prepenrichr_dlg(genelist, ...
                        backgroundlist, questtxt, parentfig)

if nargin < 4, parentfig = []; end

enrichrtype = [];
outgenelist = [];
outbackgroundlist = [];

answer = questdlg(questtxt);
if ~strcmp(answer, 'Yes'), return; end

[kgene, ybkgr, methd] = inputDialogOldVersion(parentfig);
kgene;
ybkgr;
methd;

if isempty(kgene) || isempty(ybkgr) || isempty(methd), return; end
outgenelist = genelist(1:kgene);
if ybkgr
    outbackgroundlist = backgroundlist;
end
enrichrtype = methd;

%{
answer = questdlg(sprintf('Input list contains %d genes. Run enrichment analysis with all genes?',... 
            numel(genelist)),'',...
            'Yes, use all genes', 'No, pick top k genes',...
            'Cancel', 'Yes, use all genes');
if strcmp(answer, 'Yes, use all genes')
    outgenelist = genelist;
elseif strcmp(answer, 'No, pick top k genes')
    k = gui.i_inputnumk(min([200, numel(genelist)]), 10, numel(genelist));
    if isempty(k), return; end
    outgenelist = genelist(1:k);
else
    return;
end

if ~isempty(backgroundlist)
    answer = questdlg('Enrichr with background?','');
    if strcmp(answer, 'Yes')
        outbackgroundlist = backgroundlist;
    elseif strcmp(answer, 'No')
        outbackgroundlist = []; 
    else
        return;
    end
end

enrichrtype = questdlg("Select the type of Enrichr application.","", ...
          "Web-based", "API-based", "API-based");

% answer1 = "Web-based";
% switch answer1 
%     case "API-based"
% 
% 
% 
%     case "Web-based"
%         fw = gui.gui_waitbar([], false, 'Sending genes to web browser...');
%         % gui.i_enrichtest(genelist, backgroundlist, numel(genelist));
%             if needbckg
%                 run.web_Enrichr_bkg(genelist, backgroundlist, numel(genelist));
%             else
%                 run.web_Enrichr(genelist, numel(genelist));
%             end
%         gui.gui_waitbar(fw, false, 'Check web browser & submit genes to Enrichr.');
%     otherwise
%         return;
% end
% 
%}


    function [kgene, ybkgr, methd] = inputDialogOldVersion(parentfig)

        kgene = [];
        ybkgr = [];
        methd = [];        

    if nargin < 1, parentfig = []; end
    % Create a dialog figure
    dlg = figure('Position',[500 500 300 200], 'Name', 'Input Dialog', ...
        'MenuBar', 'none', 'ToolBar', 'none', 'NumberTitle', 'off', ...
        'Resize', 'off', 'Visible', 'off');

    % Integer input field
    uicontrol(dlg, 'Style', 'text', 'Position', [10 150 120 20], 'String', 'Top k genes:');
    intInput = uicontrol(dlg, 'Style', 'edit', 'Position', [150 150 100 25], 'String', '200', 'BackgroundColor', 'white');

    % Logical checkbox input
    uicontrol(dlg, 'Style', 'text', 'Position', [10 110 120 20], 'String', 'Add background list:');
    logicalInput = uicontrol(dlg, 'Style', 'checkbox', 'Position', [150 110 100 25]);

    % Radio button group
    uicontrol(dlg, 'Style', 'text', 'Position', [10 70 120 20], 'String', 'Select Enrichr type:');
    bg = uibuttongroup('Parent', dlg, 'Position', [0.5 0.2 0.4 0.25], 'BorderType', 'none');
    rb1 = uicontrol(bg, 'Style', 'radiobutton', 'Position', [10 30 80 20], 'String', 'API-based', 'HandleVisibility', 'off');
    rb2 = uicontrol(bg, 'Style', 'radiobutton', 'Position', [10 10 80 20], 'String', 'Web-based', 'HandleVisibility', 'off');

    % OK button to confirm input
    uicontrol(dlg, 'Style', 'pushbutton', 'Position', [100 10 100 30], 'String', 'OK', 'Callback', @(src,event) processInputs());
    if ~isempty(parentfig)
        gui.i_movegui2parent(dlg, parentfig);
    end
    dlg.Visible = 'on';
    uiwait(dlg);

    % Function to process the inputs
    function processInputs()
        % Get the integer input value
        intVal = str2double(get(intInput, 'String'));
        % Get the checkbox logical value
        logicalVal = get(logicalInput, 'Value');
        % Get the selected radio button value
        selectedOption = get(get(bg, 'SelectedObject'), 'String');
        
        kgene = intVal;
        ybkgr = logicalVal==1;
        methd = selectedOption;

        % % Display the inputs
        % disp(['Integer: ', num2str(intVal)]);
        % disp(['Logical: ', num2str(logicalVal)]);
        % disp(['Selected Option: ', selectedOption]);
        
        % Close the dialog
        close(dlg);
    end
    end

end

