function [indx]=i_listdlg(a)
if nargin<1, a={'a','b','c'}; end
[indx,tf]=listdlg('ListString',a);
if tf==1
    indx;
    gui.i_listdlg(a);
elseif tf==0
    
end
end