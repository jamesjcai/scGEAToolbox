function [ok,msg]=commoncheck_R(rscriptdir)

ok=false; msg=[];

if isempty(pkg.FindRpath)
   msg=('Rscript.exe is not found');
   return;
end
%answer=questdlg('This command will check for necessary R package(s) and install them. Continue?');
%if ~strcmp(answer,'Yes')
%    msg=('Action cancelled.');
%    return;
%end
folder=fileparts(mfilename('fullpath'));
a=strfind(folder,filesep);
folder=extractBefore(folder,a(end)+1);
wrkpth=fullfile(folder,'thirdparty',rscriptdir);
cd(wrkpth);
fprintf('CURRENTWDIR = "%s"\n',wrkpth);
[~,cmdout]=pkg.RunRcode('require.R');
if strfind(cmdout,'there is no package')>0
    msg=cmdout;
    return;
end
ok=true;
end




