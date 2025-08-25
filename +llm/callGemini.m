function [done, res] = callGemini(apikey, prompt, model)
    done = false;
    % Default model if not specified
    if nargin < 3
        model = "gemini-2.0-flash";
    end
    if nargin < 1, apikey = []; end
    if nargin < 2, prompt = 'Why is the sky blue?'; end
        
    if isempty(apikey)
        preftagname = 'llapikeyenvfile';
        apikey = getpref('scgeatoolbox', preftagname);
    end

    if ~isempty(apikey) && exist(apikey,"file")
        loadenv(apikey,"FileType","env");
        apikey = getenv("GEMINI_API_KEY");
    end

    query = struct("contents",[]);
    query.contents = {struct("parts",[])};
    query.contents{1}.parts{1} = {struct("text",prompt)};    

    endpoint = "https://generativelanguage.googleapis.com/v1beta/";
    method = "generateContent";
    
    import matlab.net.*
    import matlab.net.http.*
    headers = HeaderField('Content-Type', 'application/json');
    request = RequestMessage('post', headers, query);

    fprintf('Sending request to Gemini API...\n');
    response = send(request, URI(endpoint + ...
        "models/" + model + ...
        ":" + method + ...
        "?key=" + apikey));

    if response.StatusCode == "OK"
        res = response.Body.Data.candidates.content.parts.text;
        fprintf('Response received successfully.\n');
        done = true;
    else
        res = response.Body.Data.error;
    end
    % disp(res);

    %{

    endpoint = "https://generativelanguage.googleapis.com/v1beta/";
    method = "generateContent";
    
    import matlab.net.*
    import matlab.net.http.*  
    headers = HeaderField('Content-Type', 'application/json');
    request = RequestMessage('post', headers, prompt);

    try
        response = send(request, URI(endpoint + ...
            "models/" + model + ...
            ":" + method + ...
            "?key=" + apikey));
    catch e
        fprintf('Error calling Gemini API: %s\n', e.message);
        response = [];
    end
    %}




    %{
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
    %}
end