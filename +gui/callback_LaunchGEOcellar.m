function callback_LaunchGEOcellar(src, ~)
% CALLBACK_LAUNCHGEOCELLAR  Launch the GEOcellar multi-agent pipeline UI.
%
%   Called from the scGEAToolbox main menu (Tools → GEOcellar Agent).
%   Offers three execution modes:
%     1. Completely remote  — all computation runs on the remote server
%     2. Hybrid             — Phase 1 local, Phase 2 on remote server
%     3. Completely local   — all computation runs in MATLAB

[FigureHandle, ~] = gui.gui_getfigsce(src);

if ~pkg.i_license
    gui.myErrordlg(FigureHandle, ...
        "This function requires passkey validation. You can " + ...
        "validate your passkey by selecting Help → Validate " + ...
        " Passkey from the menu.");
    return;
end

preftagname ='geocellarmodeid';
defaultindx = getpref('scgeatoolbox', preftagname, 1);

% Choose execution mode
modeOptions = [ ...
    "Remote  — server handles everything", ...
    "Hybrid  — local hypotheses, remote analysis", ...
    "Local   — full pipeline in MATLAB"];
[sel, ok] = gui.myListdlg(FigureHandle, modeOptions, ...
    "Select Execution Mode", defaultindx, false, false, [300, 180]);
if ~ok, return; end
mode = sel;   % 1 = remote, 2 = hybrid, 3 = local
setpref('scgeatoolbox', preftagname, sel);

% For hybrid and local modes the MATLAB LLM Add-On and API key are required
if mode > 1
    if ~exist("openAIChat", "file")
        gui.myHelpdlg(FigureHandle, ...
            ["The 'Large Language Models (LLMs) with MATLAB' Add-On is required." ...
             newline "Install it from the MATLAB Add-On Explorer."]);
        return;
    end

    cfg = llm.geocellar.geocellar_config();
    if isempty(cfg.OpenAIAPIKey)
        answer = gui.myQuestdlg(FigureHandle, ...
            "No API key found. Would you like to locate your llm_api_key.env file?");
        if strcmp(answer, "Yes")
            llm.i_checkllm([], [], FigureHandle);
        end
        return;
    end
end

% Set working directory — used as DataDir for downloaded GEO samples
extprogname = 'GEOcellar';
preftagname = 'externalwrkpath';
[wrkdir] = gui.gui_setprgmwkdir(extprogname, preftagname, FigureHandle);
if isempty(wrkdir), return; end

switch mode
    case 1
        llm.geocellar.geocellar_remote(wrkdir, FigureHandle);
    case 2
        llm.geocellar.geocellar_app(wrkdir, FigureHandle, UseRemote=true);
    case 3
        llm.geocellar.geocellar_app(wrkdir, FigureHandle, UseRemote=false);
end
end
