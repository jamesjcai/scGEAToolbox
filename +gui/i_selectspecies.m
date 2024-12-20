function [speciestag] = i_selectspecies(n, shorttag)

    if nargin < 2, shorttag = false; end
    if nargin < 1, n = 2; end
    
    if ~ispref('scgeatoolbox', 'preferredspecies')
        setpref('scgeatoolbox', 'preferredspecies', 'human');
    end
    preferredspecies = getpref('scgeatoolbox', 'preferredspecies', 'human');
    
    if n == 3
        answer = questdlg('Which species?', ...
            'Select Species', 'Human', 'Mouse', 'Zebrafish', preferredspecies);
    elseif n == 2
        answer = questdlg('Which species?', ...
            'Select Species', 'Human', 'Mouse', preferredspecies);
    end
    
    switch answer
        case 'Human'
            speciestag = "human";
            setpref('scgeatoolbox', 'preferredspecies', 'Human');
        case 'Mouse'
            speciestag = "mouse";
            setpref('scgeatoolbox', 'preferredspecies', 'Mouse');
        case 'Zebrafish'
            speciestag = "zebrafish";
            setpref('scgeatoolbox', 'preferredspecies', 'Zebrafish');
        otherwise
            speciestag = "";
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