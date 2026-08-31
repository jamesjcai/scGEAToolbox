function [provider, model] = i_llmsettings(overrideModel, overrideProvider)
%I_LLMSETTINGS The LLM provider and model configured for the toolbox.
%
%   [provider, model] = LLM.I_LLMSETTINGS()
%   [provider, model] = LLM.I_LLMSETTINGS(overrideModel, overrideProvider)
%
%   Reads the 'llmodelprovider' preference that GUI.I_SETLLMMODEL writes
%   (Options > Set LLM Provider & Model...), stored as "Provider:Model".
%   Either half can be overridden; pass "" to take the configured value.
%
%   The model half may itself contain colons - Ollama tags look like
%   "phi4:14b" - so only the FIRST colon separates provider from model.
%
%   Errors when nothing is configured, rather than guessing a provider the
%   user has no credentials for.
%
%   See also GUI.I_SETLLMMODEL, LLM.I_ASKLLM, LLM.I_CHECKLLM.

arguments
    overrideModel (1,1) string = ""
    overrideProvider (1,1) string = ""
end

provider = overrideProvider;
model = overrideModel;

if strlength(provider) > 0 && strlength(model) > 0
    return
end

if ~ispref('scgeatoolbox', 'llmodelprovider')
    error('llm:i_llmsettings:NotConfigured', ...
        ['No LLM provider is configured. Set one with ', ...
         'gui.i_setllmmodel, or the menu Options > Set LLM Provider & ', ...
         'Model..., then try again.']);
end

s = string(getpref('scgeatoolbox', 'llmodelprovider'));
if strlength(strtrim(s)) == 0
    error('llm:i_llmsettings:NotConfigured', ...
        'The llmodelprovider preference is empty. Re-run gui.i_setllmmodel.');
end

% Split on the FIRST colon only: "Ollama:phi4:14b" is provider "Ollama",
% model "phi4:14b".
idx = strfind(s, ":");
if isempty(idx)
    error('llm:i_llmsettings:BadPreference', ...
        ['The llmodelprovider preference is "%s", which has no ', ...
         '"Provider:Model" form. Re-run gui.i_setllmmodel.'], s);
end

if strlength(provider) == 0
    provider = extractBefore(s, idx(1));
end
if strlength(model) == 0
    model = extractAfter(s, idx(1));
end

if strlength(strtrim(model)) == 0
    error('llm:i_llmsettings:NoModel', ...
        ['Provider "%s" is configured but no model is. Re-run ', ...
         'gui.i_setllmmodel.'], provider);
end
end
