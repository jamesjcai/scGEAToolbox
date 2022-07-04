function [status]=CooccurrenceAffinity(X)
    [status]=0;
    isdebug=true;
    oldpth=pwd();
    [isok,msg]=commoncheck_R('R_CooccurrenceAffinity');
    if ~isok, error(msg); return; end

    tmpfilelist={'input.h5','output.Rds'};

    if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
    X=double(X>0);

    h5create('input.h5', '/X', size(X));
    h5write('input.h5', '/X', X);

    pkg.RunRcode('script.R');
    % [status]=copyfile('output.Rds',filename,'f');
    if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
    cd(oldpth);
end