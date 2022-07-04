function [p]=CooccurrenceAffinity(X)

isdebug=false;
oldpth=pwd();
[isok,msg]=commoncheck_R('R_CooccurrenceAffinity');
if ~isok, error(msg); return; end

tmpfilelist={'input.h5','output.h5'};

if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
X=double(X>0);

h5create('input.h5', '/X', size(X));
h5write('input.h5', '/X', X);

pkg.RunRcode('script.R');
p=h5read('output.h5','/p');
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
cd(oldpth);

end