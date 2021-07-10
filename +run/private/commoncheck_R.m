function [ok,msg]=commoncheck_R(rscriptdir)

ok=false; msg=[];

if isempty(pkg.FindRpath)
   msg=('Rscript.exe is not found');
   return;
end
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




