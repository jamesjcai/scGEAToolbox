function [y] = i_license

    y = false; % Initialize validkey to false

    preftagname ='registerpasskey';
    passkey = getpref('scgeatoolbox', preftagname, '');

    if ~isempty(passkey)
        pw1 = fileparts(mfilename('fullpath'));
        fx = 'LicenseData.mat';
        fx = fullfile(pw1, '..', 'assets', fx);
        if isfile(fx)
            S = load(fx, 'storedHash');
            storedHash = S.storedHash; 
            enteredHash = string(pkg.i_toHash(passkey));
            if enteredHash == storedHash
                y = true; % Set validkey to true if the entered hash matches the stored hash
                % disp('Passkey validated successfully.'); % Inform the user of successful validation
            end
        end
    end

%    storedHash = toHash('passkeyhere');
%    save('LicenseData.mat', 'storedHash');

end

