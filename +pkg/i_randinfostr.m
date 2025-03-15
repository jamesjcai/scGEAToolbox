function [infotagstr] = i_randinfostr(strLength)

if nargin < 1, strLength = 10; end

        chars = ['A':'Z' 'a':'z' '0':'9']; % Alphanumeric character set
        % strLength = 10; % Define desired string length
        infotagstr = chars(randperm(numel(chars), strLength));    
        
end