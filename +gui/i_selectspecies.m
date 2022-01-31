function [speciestag] = i_selectspecies

answer = questdlg('Which species?',...
    'Select Species', 'Human', 'Mouse','Zebrafish','Human');
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