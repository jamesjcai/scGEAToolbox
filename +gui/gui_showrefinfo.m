function [txt,T] = gui_showrefinfo(reftarget)

txt = [];
pw1 = fileparts(mfilename('fullpath'));
fname = fullfile(pw1, '..','resources','refinfo.txt');
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
    % uiwait(helpdlg(txt, reftarget));
    uiwait(msgbox(txt, reftarget, "help", "modal"));
end
