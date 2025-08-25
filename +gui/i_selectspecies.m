function [speciestag] = i_selectspecies(n, shorttag, parentfig)
if nargin<3, parentfig = []; end
    if nargin < 2, shorttag = false; end
    if nargin < 1, n = 2; end
    
    if ~ispref('scgeatoolbox', 'preferredspecies')
        setpref('scgeatoolbox', 'preferredspecies', 'human');
    end
    preferredspecies = getpref('scgeatoolbox', 'preferredspecies', 'human');
    
    if n == 3
        answer = gui.myQuestdlg(parentfig, 'Which species?', ...
            'Select Species',{'Human', 'Mouse', 'Zebrafish'}, preferredspecies);
    elseif n == 2
        answer = gui.myQuestdlg(parentfig, 'Which species?', ...
            'Select Species',{'Human', 'Mouse'}, preferredspecies);
    end
    speciestag = '';
    if isempty(answer), return; end
    switch answer
        case 'Human'
            speciestag = 'human';
            setpref('scgeatoolbox', 'preferredspecies', 'Human');
        case 'Mouse'
            speciestag = 'mouse';
            setpref('scgeatoolbox', 'preferredspecies', 'Mouse');
        case 'Zebrafish'
            speciestag = 'zebrafish';
            setpref('scgeatoolbox', 'preferredspecies', 'Zebrafish');
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