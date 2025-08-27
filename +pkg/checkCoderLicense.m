%% Method 1: Using license function (most reliable)
function isCoderAvailable = checkCoderLicense()
    try
        isCoderAvailable = license('test', 'MATLAB_Coder');
        if isCoderAvailable
            % fprintf('MATLAB Coder license is available.\n');
        else
            % fprintf('MATLAB Coder license is NOT available.\n');
        end
    catch
        isCoderAvailable = false;
        % fprintf('Error checking MATLAB Coder license.\n');
    end
end



%% Method 2: Using ver function to check installed toolboxes
function isCoderInstalled = checkCoderInstallation()
    toolboxes = ver;
    coderIdx = strcmpi({toolboxes.Name}, 'MATLAB Coder');
    isCoderInstalled = any(coderIdx);
    
    if isCoderInstalled
        coderInfo = toolboxes(coderIdx);
        fprintf('MATLAB Coder is installed - Version: %s\n', coderInfo.Version);
    else
        fprintf('MATLAB Coder is NOT installed.\n');
    end
end

%% Method 3: Using exist function to check for key functions
function isCoderFunctional = checkCoderFunctions()
    keyFunctions = {'codegen', 'coder.typeof', 'emlmex'};
    isCoderFunctional = true;
    
    for i = 1:length(keyFunctions)
        if exist(keyFunctions{i}, 'file') ~= 2
            fprintf('Function %s not found.\n', keyFunctions{i});
            isCoderFunctional = false;
        end
    end
    
    if isCoderFunctional
        fprintf('MATLAB Coder key functions are available.\n');
    else
        fprintf('MATLAB Coder key functions are NOT available.\n');
    end
end

%% Method 4: Comprehensive check function
function [hasLicense, isInstalled, isFunctional] = checkMATLABCoder()
    fprintf('=== MATLAB Coder Availability Check ===\n\n');
    
    % Check license
    hasLicense = checkCoderLicense();
    
    % Check installation
    isInstalled = checkCoderInstallation();
    
    % Check functionality
    isFunctional = checkCoderFunctions();
    
    % Summary
    fprintf('\n=== Summary ===\n');
    fprintf('License Available: %s\n', logical2str(hasLicense));
    fprintf('Toolbox Installed: %s\n', logical2str(isInstalled));
    fprintf('Functions Available: %s\n', logical2str(isFunctional));
    
    if hasLicense && isInstalled && isFunctional
        fprintf('✓ MATLAB Coder is fully available and ready to use.\n');
    else
        fprintf('✗ MATLAB Coder has some limitations or is not available.\n');
    end
end

%% Helper function to convert logical to string
function str = logical2str(val)
    if val
        str = 'Yes';
    else
        str = 'No';
    end
end

%% Quick one-liner check
function result = isCoderAvailable()
    result = license('test', 'MATLAB_Coder') && exist('codegen', 'file') == 2;
end

%% Usage examples:
% checkMATLABCoder()  % Comprehensive check
% isCoderAvailable()  % Quick boolean check