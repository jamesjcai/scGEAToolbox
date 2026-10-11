function [done, res] = callClaude(prompt, model)
%CALLCLAUDE Send one prompt to the Anthropic Messages API (Claude).
%   [done, res] = llm.callClaude(prompt, model)
%
%   prompt  text to send (default: 'Why is the sky blue?')
%   model   model name (default: "claude-sonnet-5-5")
%
%   done is true when the API answered OK; res is the reply text, or the
%   API's error struct (with a .message field) when it did not.
%
%   The key is ANTHROPIC_API_KEY, looked up by llm.i_getapikey: the
%   environment, the MATLAB vault, then the env file named by the
%   'llapikeyenvfile' preference. Store it once with
%   llm.i_storeapikey("ANTHROPIC_API_KEY").
%
%   Talks to api.anthropic.com directly. Claude is also reachable through
%   the TAMU proxy (llm.callOpenAIChatTAMU with a "protected.Claude ..."
%   model), which needs no Anthropic key.
%
%   See also llm.i_askllm, llm.callGemini, llm.callOpenAIChat.
done = false;
if nargin < 2, model = "claude-sonnet-5-5"; end
if nargin < 1, prompt = 'Why is the sky blue?'; end

apikey = llm.i_getapikey("ANTHROPIC_API_KEY", llm.i_readkeyfile());
if isempty(apikey)
    error('llm:callClaude:noKey', ...
        ['ANTHROPIC_API_KEY is not set. Store it with ' ...
        'llm.i_storeapikey("ANTHROPIC_API_KEY"), or put it in the env file ' ...
        'named by the ''llapikeyenvfile'' preference.']);
end

% Long enough for the LLM report summaries; a label answer stops far short.
maxTokens = 8192;
query = struct("model", string(model), "max_tokens", maxTokens, ...
    "messages", {{struct("role", "user", "content", char(prompt))}});

import matlab.net.*
import matlab.net.http.*
headers = [HeaderField('Content-Type', 'application/json'), ...
    HeaderField('x-api-key', apikey), ...
    HeaderField('anthropic-version', '2023-06-01')];
request = RequestMessage('post', headers, query);
% HTTPOptions' ResponseTimeout is Inf by default, which would let a stalled
% request hang the caller.
opts = HTTPOptions('ConnectTimeout', 30, 'ResponseTimeout', 300);

fprintf('Sending request to Anthropic API...\n');
response = send(request, URI("https://api.anthropic.com/v1/messages"), opts);
data = response.Body.Data;

if response.StatusCode == "OK"
    % Every text block of the reply. JSONDECODE gives a cell array when the
    % blocks differ in fields (e.g. a thinking block before the text).
    blocks = data.content;
    if ~iscell(blocks)
        blocks = num2cell(blocks);
    end
    isText = cellfun(@(b) isfield(b, 'type') && strcmp(b.type, 'text'), blocks);
    texts = cellfun(@(b) string(b.text), blocks(isText), 'UniformOutput', false);
    res = strjoin(string([texts{:}]), "");
    if isfield(data, 'stop_reason') && strcmp(data.stop_reason, 'max_tokens')
        fprintf('Note: the reply was cut off at %d tokens.\n', maxTokens);
    end
    fprintf('Response received successfully.\n');
    done = true;
elseif isstruct(data) && isfield(data, 'error')
    res = data.error;
    disp(res);
else
    % A proxy or gateway error page has no JSON error body.
    res = struct("type", "http_error", "message", ...
        "HTTP " + string(response.StatusCode) + " from api.anthropic.com");
    disp(res);
end
end
