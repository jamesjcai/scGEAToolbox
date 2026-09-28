function [done, res] = callOpenAIChatTAMU(prompt, model)
%CALLOPENAICHATTAMU One prompt to the TAMU AI chat API.
%   [done, res] = llm.callOpenAIChatTAMU(prompt, model) reads TAMUAI_API_KEY
%   and TAMUAI_API_BASE (default https://chat-api.tamu.ai/api).
%   See also llm.callOpenAIChat, llm.i_askllm.
if nargin < 2, model = "protected.gemini-2.0-flash-lite"; end
if nargin < 1, prompt = 'Why is the sky blue?'; end
[done, res] = i_callopenaicompat(prompt, model, "TAMUAI_API_KEY", ...
    "TAMUAI_API_BASE", "https://chat-api.tamu.ai/api", "TAMU AI");
end
