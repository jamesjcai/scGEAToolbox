function response = callGemini(apiKey, prompt, model)
    % Default model if not specified
    if nargin < 3, model = 'gemini-2.0-flash'; end
    if nargin < 1, apiKey = []; end
    if nargin < 2, prompt = 'Why is the sky blue?'; end

    if isempty(apiKey)
        preftagname = 'llapikeyenvfile';
        apiKey = getpref('scgeatoolbox', preftagname);
    end

    if ~isempty(apiKey) && exist(apiKey,"file")
        loadenv(apiKey,"FileType","env");
        apiKey = getenv("GEMINI_API_KEY");
    end
    
    % API endpoint
    url = ['https://generativelanguage.googleapis.com/v1beta/models/' model ':generateContent?key=' apiKey];
          % https://generativelanguage.googleapis.com/v1beta/models/gemini-2.0-flash:generateContent?key=
           % https://generativelanguage.googleapis.com/v1beta/models/gemini-2.0-flash:generateContent?key=AIzaSyDs2iibCA_5eeEgixAH8tzjjDK1j_CM_Nk
    % Construct the request body
    requestBody = struct('contents', struct('parts', struct('text', prompt)));
    jsonBody = jsonencode(requestBody);
    
    % Set up the HTTP request
    options = weboptions('ContentType', 'json', 'RequestMethod', 'post');
    
    % Make the API call
    try
        response = webwrite(url, jsonBody, options);
        % Extract the actual response text
        if isfield(response, 'candidates') && ~isempty(response.candidates) && ...
           isfield(response.candidates{1}, 'content') && ...
           isfield(response.candidates{1}.content, 'parts') && ...
           ~isempty(response.candidates{1}.content.parts)
            response = response.candidates{1}.content.parts{1}.text;
        end
    catch e
        fprintf('Error calling Gemini API: %s\n', e.message);
        response = [];
    end
end