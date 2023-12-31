function [answer]=i_questdlg(quest,btn1,btn2,btn3,defbtn)
if nargin < 5, defbtn = btn1; end
answer = questdlg(quest,'',btn1,btn2,btn3,defbtn);
switch answer
    case btn1
        answer = 1;
    case btn2
        answer = 2;
    case btn3
        answer = 3;
    otherwise
        answer = [];
end