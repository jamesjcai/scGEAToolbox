function [done, res] = callGemini(apikeyfile, prompt, model)
%CALLGEMINI Send one prompt to the Gemini generateContent API.
%   [done, res] = llm.callGemini(apikeyfile, prompt, model)
%
%   apikeyfile  env file holding GEMINI_API_KEY (default: the scgeatoolbox
%               'llapikeyenvfile' preference). If it is empty or not a file,
%               GEMINI_API_KEY is read from the environment as it stands.
%   prompt      text to send (default: 'Why is the sky blue?')
%   model       model name (default: "gemini-2.5-flash")
%
%   done is true when the API answered OK; res is the reply text, or the
%   API's error struct when it did not.
%
%   The key travels in the x-goog-api-key header. It used to be appended to
%   the URL as ?key=..., which puts it in proxy and server access logs.
done = false;
% Default model if not specified
if nargin < 3
    model = "gemini-2.5-flash";
end
if nargin < 1, apikeyfile = []; end
if nargin < 2, prompt = 'Why is the sky blue?'; end

if isempty(apikeyfile) && ispref('scgeatoolbox', 'llapikeyenvfile')
    apikeyfile = getpref('scgeatoolbox', 'llapikeyenvfile');
end

if ~isempty(apikeyfile) && isfile(apikeyfile)
    loadenv(apikeyfile, "FileType", "env");
end
apikey = getenv("GEMINI_API_KEY");
if isempty(apikey)
    error('llm:callGemini:noKey', ...
        ['GEMINI_API_KEY is not set. Put it in the env file named by the ' ...
        '''llapikeyenvfile'' preference, or pass that file as APIKEYFILE.']);
end

% "parts" is an array of part objects, [{"text": ...}]. This used to put a
% cell inside the cell, which encodes as [[{...}]].
query = struct("contents", {{struct("parts", {{struct("text", char(prompt))}})}});

endpoint = "https://generativelanguage.googleapis.com/v1beta/";
method = "generateContent";

import matlab.net.*
import matlab.net.http.*
headers = [HeaderField('Content-Type', 'application/json'), ...
    HeaderField('x-goog-api-key', apikey)];
request = RequestMessage('post', headers, query);
% HTTPOptions' ResponseTimeout is Inf by default: a stalled request hung
% the caller, and with it the app.
opts = HTTPOptions('ConnectTimeout', 30, 'ResponseTimeout', 300);

fprintf('Sending request to Gemini API...\n');
response = send(request, URI(endpoint + "models/" + model + ":" + method), opts);

if response.StatusCode == "OK"
    % Every text part of the first candidate, not only the first part.
    parts = response.Body.Data.candidates(1).content.parts;
    if iscell(parts)
        parts = [parts{:}];
    end
    res = strjoin(string({parts.text}), "");
    fprintf('Response received successfully.\n');
    done = true;
else
    res = response.Body.Data.error;
    disp(res);
end
end
