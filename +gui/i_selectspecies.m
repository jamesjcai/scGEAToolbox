function [speciestag] = i_selectspecies(n, shorttag)

if nargin < 2, shorttag = false; end
if nargin < 1, n = 2; end

if n == 3
    answer = questdlg('Which species?', ...
        'Select Species', 'Human', 'Mouse', 'Zebrafish', 'Human');
elseif n == 2
    answer = questdlg('Which species?', ...
        'Select Species', 'Human', 'Mouse', 'Human');
end

switch answer
    case 'Human'
        speciestag = "human";
    case 'Mouse'
        speciestag = "mouse";
    case 'Zebrafish'
        speciestag = "zebrafish";
    otherwise
        speciestag = [];
        helpdlg('Action cancelled.','');
        return;
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