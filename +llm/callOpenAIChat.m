function [done, res] = callOpenAIChat(prompt, model)
%CALLOPENAICHAT One prompt to the OpenAI chat API.
%   [done, res] = llm.callOpenAIChat(prompt, model) reads OPENAI_API_KEY (or
%   OpenAI_API_KEY) from the 'llapikeyenvfile' env file or the environment.
%   See also llm.callOpenAIChatTAMU, llm.callOpenAIChatNVIDIA, llm.i_askllm.
if nargin < 2, model = "gpt-4.1"; end
if nargin < 1, prompt = 'Why is the sky blue?'; end
[done, res] = i_callopenaicompat(prompt, model, ["OPENAI_API_KEY", "OpenAI_API_KEY"], ...
    "", "", "OpenAI");
end
