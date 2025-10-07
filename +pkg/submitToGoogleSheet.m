function status = submitToGoogleSheet(username, hostname, os, matlab_version, app_version, location, comment)
% submitToGoogleSheet - Submit user information to Google Spreadsheet web app
%
% Syntax:
%   status = submitToGoogleSheet(username, hostname, os, matlab_version, app_version, location, comment)
%
% Inputs:
%   username       - String: User's username
%   hostname       - String: Computer hostname
%   os             - String: Operating system (e.g., 'Windows 10', 'macOS', 'Linux')
%   matlab_version - String: MATLAB version (e.g., 'R2023b')
%   app_version    - String: Application version
%   location       - String: User's location
%   comment        - String: Optional comment (can be empty string '')
%
% Output:
%   status - Struct with fields:
%            .success - Boolean indicating if submission was successful
%            .message - String with response message
%
% Example:
%   status = pkg.submitToGoogleSheet('john_doe', 'DESKTOP-ABC123', 'Windows 10', ...
%                                'R2023b', '1.0.0', 'New York', 'Test submission');
%
%   if status.success
%       disp('Data submitted successfully!');
%   else
%       disp(['Error: ' status.message]);
%   end

    % ========================================================================
    % CONFIGURATION: Replace with your Google Apps Script web app URL
    % ========================================================================
    WEB_APP_URL = 'https://script.google.com/macros/s/AKfycbyx-1QGE8p9Ghd6OktkeomTGaODFlQfUefL-RAC--S-7XmkroJ_OQugfzJZ4o2mqMfN_Q/exec';
    % Example: 'https://script.google.com/macros/s/AKfycby.../exec'
    
    % Validate inputs
    if nargin < 7
        error('All 7 arguments are required: username, hostname, os, matlab_version, app_version, location, comment');
    end
    
    % Convert all inputs to strings if they aren't already
    username = string(username);
    hostname = string(hostname);
    os = string(os);
    matlab_version = string(matlab_version);
    app_version = string(app_version);
    location = string(location);
    comment = string(comment);
    
    % Create data structure
    data = struct(...
        'username', char(username), ...
        'hostname', char(hostname), ...
        'os', char(os), ...
        'matlab_version', char(matlab_version), ...
        'app_version', char(app_version), ...
        'location', char(location), ...
        'comment', char(comment) ...
    );
    
    % Convert to JSON
    jsonData = jsonencode(data);
    
    % Set up web options
    options = weboptions(...
        'MediaType', 'application/json', ...
        'RequestMethod', 'post', ...
        'Timeout', 30, ...
        'ContentType', 'json' ...
    );
    
    % Initialize status structure
    status = struct('success', false, 'message', '');
    
    try
        % Send POST request
        response = webwrite(WEB_APP_URL, jsonData, options);
        
        % Parse response
        if isstruct(response)
            if isfield(response, 'status') && strcmp(response.status, 'success')
                status.success = true;
                status.message = response.message;
            else
                status.success = false;
                status.message = response.message;
            end
        else
            status.success = true;
            status.message = 'Data submitted successfully';
        end
        
    catch ME
        % Handle errors
        status.success = false;
        status.message = sprintf('Error: %s', ME.message);
    end
end