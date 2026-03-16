function [speciestag] = i_selectspecies(n, shorttag, parentfig, preferredspecies)
    if nargin < 3, parentfig = []; end
    if nargin < 2, shorttag = false; end
    if nargin < 1, n = 2; end
    
    if ~ispref('scgeatoolbox', 'preferredspecies')
        setpref('scgeatoolbox', 'preferredspecies', 'human');
    end
    if nargin < 4 || isempty(preferredspecies)
        preferredspecies = getpref('scgeatoolbox', 'preferredspecies', 'human');
    end

    if n == 3
        answer = gui.myQuestdlg(parentfig, 'Which species?', ...
            'Select Species',{'Human', 'Mouse', 'Zebrafish'}, capitalizeFirst(preferredspecies));
    elseif n == 2
        answer = gui.myQuestdlg(parentfig, 'Which species?', ...
            'Select Species',{'Human', 'Mouse'}, capitalizeFirst(preferredspecies));
    end
    speciestag = '';
    if isempty(answer), return; end
    switch answer
        case 'Human'
            speciestag = 'human';
            setpref('scgeatoolbox', 'preferredspecies', 'human');
        case 'Mouse'
            speciestag = 'mouse';
            setpref('scgeatoolbox', 'preferredspecies', 'mouse');
        case 'Zebrafish'
            speciestag = 'zebrafish';
            setpref('scgeatoolbox', 'preferredspecies', 'zebrafish');
        otherwise
            speciestag = '';
    end
    
    if shorttag
        switch lower(speciestag)
            case 'human'
                speciestag = 'hs';
            case 'mouse'
                speciestag = 'mm';
        end
    end

end


function str = capitalizeFirst(str)
    if ischar(str)
        if ~isempty(str)
            str(1) = upper(str(1));
        end
    elseif isstring(str)
        mask = strlength(str) > 0;
        str(mask) = extractBefore(str(mask), 2).upper() + extractAfter(str(mask), 1);
    end
end
