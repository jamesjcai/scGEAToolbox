function [txt, ok] = i_askllm(prompt, provider, model, quiet)
%I_ASKLLM Send one prompt to the configured LLM, whichever provider that is.
%
%   [txt, ok] = LLM.I_ASKLLM(prompt)
%   [txt, ok] = LLM.I_ASKLLM(prompt, provider, model)
%   [txt, ok] = LLM.I_ASKLLM(prompt, provider, model, quiet)
%
%   With no provider/model, both come from LLM.I_LLMSETTINGS, i.e. from what
%   GUI.I_SETLLMMODEL configured. OK is false when the call failed, in which
%   case TXT carries the reason rather than an answer - callers that retry
%   should test OK rather than parsing TXT.
%
%   Reasoning traces (<think>...</think>) are stripped from the reply, since
%   no caller wants them and every one of them used to strip them inline.
%
%   QUIET (default true) suppresses the per-call chatter the underlying
%   clients print. A caller asking once wants to see it; one asking a hundred
%   times in a loop does not.
%
%   Each LLM.CALL* client reads its own credentials from the
%   'llapikeyenvfile' preference, so nothing needs passing here. Providers
%   that GUI.I_SETLLMMODEL offers but that have no client yet raise an
%   identified error rather than failing obscurely.
%
%   See also LLM.I_LLMSETTINGS, LLM.CALLOLLAMA, LLM.CALLOPENAICHAT,
%   LLM.CALLOPENAICHATTAMU, LLM.CALLOPENAICHATNVIDIA, LLM.CALLGEMINI.

arguments
    prompt (1,1) string
    provider (1,1) string = ""
    model (1,1) string = ""
    quiet (1,1) logical = true
end

if strlength(provider) == 0 || strlength(model) == 0
    [provider, model] = llm.i_llmsettings(model, provider);
end

% ok is the guard: it stays false on the early-return paths below.
ok = false;

try
    if quiet
        % The clients print progress on every call; capture it so a loop of
        % a hundred prompts does not bury the caller's own output.
        [~, done, res] = evalc('i_dispatch(prompt, provider, model)');
    else
        [done, res] = i_dispatch(prompt, provider, model);
    end
catch ME
    txt = string(ME.message);
    return
end

if ~done
    txt = "The " + provider + " client reported a failed request.";
    return
end

% Reasoning models wrap their working in <think>...</think>. No caller in
% the toolbox wants that in a summary or a label, and six of them stripped it
% inline with their own copy of this regexp before this function existed.
txt = regexprep(string(res), '<think>[\s\S]*?</think>', '');
txt = strtrim(txt);
ok = strlength(txt) > 0;
if ~ok
    txt = "The " + provider + " client returned an empty response.";
end
end


function [done, res] = i_dispatch(prompt, provider, model)
% One arm per provider that has a working client. The clients are already
% uniform - [done, res] = call*(prompt, model) - except Gemini, which takes
% the key file first.
switch lower(provider)
    case "ollama"
        [done, res] = llm.callOllama(prompt, model);
    case "openai"
        [done, res] = llm.callOpenAIChat(prompt, model);
    case "tamuaichat"
        [done, res] = llm.callOpenAIChatTAMU(prompt, model);
    case "nvidia"
        [done, res] = llm.callOpenAIChatNVIDIA(prompt, model);
    case "gemini"
        keyfile = "";
        if ispref('scgeatoolbox', 'llapikeyenvfile')
            keyfile = string(getpref('scgeatoolbox', 'llapikeyenvfile'));
        end
        [done, res] = llm.callGemini(keyfile, prompt, model);
    otherwise
        error('llm:i_askllm:UnsupportedProvider', ...
            ['Provider "%s" has no client in +llm yet. Supported: ', ...
             'Ollama, OpenAI, TAMUAIChat, NVIDIA, Gemini. Choose one of ', ...
             'those with gui.i_setllmmodel, or pass Provider= explicitly.'], ...
            provider);
end
end
