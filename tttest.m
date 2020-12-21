pw1=fileparts(which(mfilename))
pw2=fileparts(mfilename('fullpath'))
mfilename()

pth1=fullfile(pw1,'thirdparty','R_MAST')
pth2=fullfile(pw2,'thirdparty/R_MAST')
pth1=fullfile(pw1,'../','thirdparty','R_MAST')


a=strfind(pw1,filesep)
extractBefore(pw1,a(2)+1);
