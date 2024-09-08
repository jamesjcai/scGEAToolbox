function [outgenelist, outbackgroundlist, enrichrtype] = gui_prepenrichr_dlg(genelist, backgroundlist, questtxt)

enrichrtype = [];
outgenelist = [];
outbackgroundlist = [];

answer = questdlg(questtxt);
if ~strcmp(answer, 'Yes'), return; end

    definput = {'0.05', '1.0', '0.01'};
    prompt = {'Min. abs(diff(pct)):', ...
              'Min. abs(log2(FC)):', ...
              'Adjusted P-value cutoff:'};
    dlgtitle = 'DE Result Filter';
    dims = [1, 80];
    answer = inputdlg(prompt, dlgtitle, dims, definput);

    if isempty(answer), return; end
    try
        mindiffpct = str2double(answer{1});
        minabsolfc = str2double(answer{2});
        apvaluecut = str2double(answer{3});
        assert((mindiffpct > 0) && (mindiffpct < 1));
        assert((minabsolfc >= 1) && (minabsolfc < 100));
        assert((apvaluecut >= 0.0) && (apvaluecut <= 1.0));
    catch
        errordlg('Invalid input.');
        return;
    end






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
