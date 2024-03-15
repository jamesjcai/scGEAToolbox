function i_addwd2script(scriptfile, targetdir, ftype)

    if nargin < 3, ftype = ''; end   % 'R' or 'python'
    if ~exist(scriptfile, 'file'), error('pkg.i_addwd2script error.'); end
    if ~exist(targetdir, 'dir'), error('pkg.i_addwd2script error.'); end
    
    S = fileread(scriptfile);
    
    switch ftype
        case 'R'
            x = sprintf('setwd("%s")', strrep(targetdir,'\','\\'));
            S = [x, newline, S];
            outscriptfile = fullfile(targetdir,'script.R');
        case 'python'
            x = sprintf('os.chdir("%s")', strrep(targetdir,'\','\\'));
            S = ['import os', newline, x, newline, S];
            outscriptfile = fullfile(targetdir,'script.py');
    end
    
    
    FID = fopen(outscriptfile, 'w');
    if FID == -1, error('Cannot open file %s', outscriptfile); end
    fwrite(FID, S, 'char');
    fclose(FID);
    
    
    % https://www.mathworks.com/matlabcentral/answers/80529-append-data-at-the-beginning-of-a-file-matlab
end
