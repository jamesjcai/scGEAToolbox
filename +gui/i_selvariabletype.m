function [answer] = i_selvariabletype(y,parentfig)
if nargin<2, parentfig = []; end

[c] = grp2idx(y);
n = max(c);
if n < 20
    deft = 'Categorical/Discrete';
else
    deft = 'Numerical/Continuous';
end
answer = gui.myQuestdlg(parentfig, 'What is the variable type?', '', ...
    {'Categorical/Discrete', ...
    'Numerical/Continuous', 'Unknown'}, deft);


% n=max(c);
% if n<40
%     f=0.5*(n-1)./n;
%     f=1+f.*(1:2:2*n);
%     cb=colorbar('Ticks',f,'TickLabels',cellstr(cL));
% else
%     %c=thisc;
%     set(h,'CData',thisc);
%     cb=colorbar;
%     %cb=colorbar('Ticks',[]);
% end
