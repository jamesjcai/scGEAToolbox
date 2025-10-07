function geatoolUsageLogger(src, evt, app_version)
    %=============================
    % scGEATool Usage Logger
    %=============================
    
    if nargin < 1, src = []; end
    if nargin < 2, evt = []; end
    if nargin < 3, app_version = '1.0.0'; end
    
    %% --- Metadata
    username = getenv('USERNAME');
    if isempty(username)
        username = getenv('USER');
    end
    hostname = char(java.net.InetAddress.getLocalHost.getHostName());
    os = system_dependent('getos');
    matlab_version = version;
    
    comment = ''; % optional notes
    
    %% --- Optional: IP-based Location
    location = 'Unknown';
    try
        locData = webread('https://ipapi.co/json/');
        if isfield(locData,'city') && isfield(locData,'country_name')
            location = sprintf('%s, %s', locData.city, locData.country_name);
        end
    catch
        % keep location as 'Unknown' if API fails
    end
    
%   status = pkg.submitToGoogleSheet('john_doe', 'DESKTOP-ABC123', 'Windows 10', ...
%                                'R2023b', '1.0.0', 'New York', 'Test submission');

   status = pkg.submitToGoogleSheet(username, hostname, os, ...
                                matlab_version, app_version, location, comment);

%{
    %% --- Data to Send
    data = struct(...
        'username', username, ...
        'hostname', hostname, ...
        'os', os, ...
        'matlab_version', matlab_version, ...
        'app_version', app_version, ...
        'location', location, ...
        'comment', comment ...
    );
    
    %% --- Google Apps Script URL
    % IMPORTANT: Use the deployed Web App /exec URL, NOT the "preview" URL
    % Example:
    scriptURL = 'https://script.google.com/macros/s/AKfycbwZ3FrptB2Pxr3SO4H1J_nnDT-1FdpYD1x9lP5RUOeZcD8IE1Va3aSp4COm4Rku0zpPmg/exec'; 
    
    %% --- Upload Usage Data
    try
        options = weboptions('MediaType','application/json','RequestMethod','post');
        response = webwrite(scriptURL, data, options);
        fprintf('Usage info uploaded successfully to Google Sheet.\n');
    catch ME
        % warning('Could not sync with Google Sheet: %s', ME.message);
    end    
%}
    
end
