function [done, res] = callOpenAIChatNVIDIA(prompt, model)
%CALLOPENAICHATNVIDIA One prompt to the NVIDIA chat API.
%   [done, res] = llm.callOpenAIChatNVIDIA(prompt, model) reads
%   NVIDIA_API_KEY and NVIDIA_API_BASE (default
%   https://integrate.api.nvidia.com/v1).
%   See also llm.callOpenAIChat, llm.i_askllm.
if nargin < 2, model = "minimaxai/minimax-m2.5"; end
if nargin < 1, prompt = 'Why is the sky blue?'; end
[done, res] = i_callopenaicompat(prompt, model, "NVIDIA_API_KEY", ...
    "NVIDIA_API_BASE", "https://integrate.api.nvidia.com/v1", "NVIDIA");
end
