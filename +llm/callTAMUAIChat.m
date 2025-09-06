function [done, response] = callTAMUAIChat(apikey, prompt, model)
    
    done = false;
    % CALLTAMUCHAT Call TAMU Chat API for text generation
    %
    % Syntax:
    %   response = callTAMUChat(apiKey, prompt)
    %   response = callTAMUChat(apiKey, prompt, model)
    %
    % Inputs:
    %   apiKey - String, API key for authentication
    %   prompt - String, the user prompt/question
    %   model  - String, model name (optional, default: 'protected.llama3.2')
    %
    % Output:
    %   response - String, the generated response text, or empty if error
    %
    % Example:
    %   apiKey = "";
    %   response = callTAMUChat(apiKey, "Why is the sky blue?");
    
    % Default model if not specified
    if nargin < 3, model = 'protected.gpt-4.1'; end
    if nargin < 1, apikey = []; end
    if nargin < 2, prompt = 'Why is the sky blue?'; end

    if isempty(apikey)
        preftagname = 'llapikeyenvfile';
        apikey = getpref('scgeatoolbox', preftagname);
    end

    if ~isempty(apikey) && exist(apikey,"file")
        loadenv(apikey,"FileType","env");
        apikey = getenv("TAMUAI_API_KEY");
    end    
    
    % API endpoint
    url = 'https://chat-api.tamu.ai/api/chat/completions';
    
    % Construct the request body
    messages_cell = {struct('role', 'user', 'content', prompt)};
    requestBody = struct('model', model, ...
                        'stream', false, ...
                        'messages', {messages_cell});
    
    % Set up the HTTP request options
    options = weboptions('MediaType', 'application/json', ...
                        'RequestMethod', 'POST', ...
                        'Timeout', 30, ...
                        'HeaderFields', {'Authorization', sprintf('Bearer %s', apikey)});
    
    % Make the API call
    try
        fprintf('Sending request to TAMU Chat API...\n');
        apiResponse = webwrite(url, requestBody, options);
        
        % Extract the actual response text
        if isfield(apiResponse, 'choices') && ~isempty(apiResponse.choices)
            % Handle both cell array and struct array cases
            if iscell(apiResponse.choices)
                choice = apiResponse.choices{1};
            else
                choice = apiResponse.choices(1);
            end
            
            if isfield(choice, 'message') && isfield(choice.message, 'content')
                response = choice.message.content;
                fprintf('Response received successfully.\n');
                done = true;
            else
                fprintf('Unexpected choice structure in TAMU Chat API response\n');
                fprintf('Choice structure:\n');
                disp(choice);
                response = [];
            end
        else
            fprintf('No choices field found in TAMU Chat API response\n');
            fprintf('Full response structure:\n');
            disp(apiResponse);
            response = [];
        end
        
    catch e
        if contains(e.message, 'timed out')
            fprintf('Error: Request timed out. The TAMU Chat API server may be slow or unavailable.\n');
            fprintf('You can try:\n');
            fprintf('  1. Running the request again\n');
            fprintf('  2. Checking if the API endpoint is accessible\n');
            fprintf('  3. Increasing the timeout further if needed\n');
        else
            fprintf('Error calling TAMU Chat API: %s\n', e.message);
        end
        response = [];
    end
end