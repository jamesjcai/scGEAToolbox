function callback_Passkey(src, ~)

    [FigureHandle, ~] = gui.gui_getfigsce(src);
    if pkg.i_license
        gui.myHelpdlg(FigureHandle, 'Passkey validated successfully.');    
    else
        done = false;
        if gui.i_isuifig(FigureHandle)
            % newkey = gui.myInputdlg({'Passkey'}, 'Input', {''}, FigureHandle);
            newkey = gui.i_inputdlg('Passkey','', FigureHandle);
        else
            newkey = inputdlg('Passkey', 'Input', [1, 50], {''});
        end
        
        if isempty(newkey)
            return;
        else
            newkey = cell2mat(newkey);
            if ~isempty(newkey)
                pw1 = fileparts(mfilename('fullpath'));
                fx = 'LicenseData.mat';
                fx = fullfile(pw1, '..', 'assets', fx);
                if isfile(fx)
                    S = load(fx, 'storedHash');
                    
                    storedHash = string(S.storedHash);
                    enteredHash = string(pkg.i_toHash(newkey));

                    if enteredHash == storedHash
                        preftagname ='registerpasskey';
                        setpref('scgeatoolbox', preftagname, newkey);
                        gui.myHelpdlg(FigureHandle, 'Passkey validated successfully.');
                        done = true;
                    end
                end
            end
        end
        if ~done
            gui.myErrordlg(FigureHandle, 'Passkey is not validated.');
        end
    end

end