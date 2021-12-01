function [speciestag] = i_selectspecies

answer = questdlg('Which species?',...
    'Select Species', 'Mouse', 'Human', 'Mouse');
switch answer
    case 'Human'
        speciestag = "human";
    case 'Mouse'
        speciestag = "mouse";
    otherwise
        speciestag = [];
end
end