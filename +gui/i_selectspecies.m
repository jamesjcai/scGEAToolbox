function [speciestag] = i_selectspecies(n)

if nargin<1, n=3; end
if n==3
    answer = questdlg('Which species?',...
        'Select Species', 'Human', 'Mouse','Zebrafish','Human');
elseif n==2
    answer = questdlg('Which species?',...
        'Select Species', 'Human', 'Mouse','Human');
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
end
end