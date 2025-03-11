function [y, txt, T] = gui_showrefinfo(reftarget, FigureHandle)

%see also: gui.gui_uishowrefinfo
if nargin<2, FigureHandle = []; end
y=false;
txt = [];
pw1 = fileparts(mfilename('fullpath'));
fname = fullfile(pw1, '..','resources','Misc','refinfo.txt');
fid=fopen(fname,'r');
T=textscan(fid,'%s%s','Delimiter','\t');
fclose(fid);
reftag=string(T{:,1});

idx=find(reftarget==reftag);
if ~isempty(idx)
    txt=T{:,2}{idx};
end

if ~isempty(txt)
    fprintf('%s\n%s\n', reftarget, txt);
    answer = gui.myQuestdlg(FigureHandle, ...
         txt, reftarget, {'Continue','Cancel'}, 'Continue');
     switch answer
         case 'Continue'
             y = true;
         otherwise
     end
end
