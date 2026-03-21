function parameters = buildOpenAIParameters(messages, functions, nvp)
%buildOpenAIParameters Build parameters struct for OpenAI Chat API
%
%   PARAMETERS = buildOpenAIParameters(MESSAGES, FUNCTIONS, NVP) builds a
%   struct in the format expected by the OpenAI Chat Completions API,
%   combining MESSAGES, FUNCTIONS and parameters in NVP.
%
%   NVP is a struct with fields: ToolChoice, ModelName, Temperature, TopP,
%   NumCompletions, StopSequences, MaxNumTokens, PresencePenalty,
%   FrequencyPenalty, ResponseFormat, Seed, StreamFun.
%
%   See also: llms.internal.callOpenAIChatAPI

%   Copyright 2026 The MathWorks, Inc.

parameters = llms.internal.buildAzureParameters(messages, functions, nvp);

parameters.model = nvp.ModelName;

end
