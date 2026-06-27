function callback_LaunchGEOcellar(src, ~)
% CALLBACK_LAUNCHGEOCELLAR  Launch the GEOcellar multi-agent pipeline UI.
%
%   Called from the scGEAToolbox main menu (Tools → GEOcellar Agent).
%   Offers four execution modes:
%     1. Local        — full pipeline in MATLAB
%     2. Hybrid (P1)  — local hypotheses, remote analysis
%     3. Hybrid (P2)  — remote hypotheses, local analysis
%     4. Remote       — server handles everything

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
    "Local          — full pipeline in MATLAB", ...
    "Hybrid (P1)    — local hypotheses, remote analysis", ...
    "Hybrid (P2)    — remote hypotheses, local analysis", ...
    "Remote         — server handles everything"];
[sel, ok] = gui.myListdlg(FigureHandle, modeOptions, ...
    "Select Execution Mode", defaultindx, false, false, [340, 210]);
if ~ok, return; end
mode = sel;   % 1 = local, 2 = P1 local/P2 remote, 3 = P1 remote/P2 local, 4 = remote
setpref('scgeatoolbox', preftagname, sel);

% Modes 1 and 2 run Phase 1 locally — require the MATLAB LLM Add-On and API key
if mode <= 2
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
        llm.geocellar.geocellar_app(wrkdir, FigureHandle);
    case 2
        llm.geocellar.geocellar_app(wrkdir, FigureHandle, UseRemote=true);
    case 3
        llm.geocellar.geocellar_app(wrkdir, FigureHandle, UseRemote1=true);
    case 4
        llm.geocellar.geocellar_remote(wrkdir, FigureHandle);
end
end
