function [ok,wrkpth,x]=pycommon(prgwkdir)
arguments
    prgwkdir {mustBeTextScalar}
end

    ok=false;

    oldpth=pwd();
    pw1=fileparts(mfilename('fullpath'));
    wrkpth=fullfile(pw1,'external',prgwkdir);
    cd(wrkpth);    
       
    x=pyenv;
    if strlength(x.Executable)==0, return; end

    fw = gui.gui_waitbar([],[],'Checking Python environment...');

    try
        pkg.i_add_conda_python_path;
    catch
        
    end
    cmdlinestr=sprintf('"%s" "%s%srequire.py"', ...
            x.Executable,wrkpth,filesep);
    disp(cmdlinestr)
    [status,cmdout]=system(cmdlinestr,'-echo');
    if status~=0
        cd(oldpth);
        if isvalid(fw) 
            gui.gui_waitbar(fw,[],'Checking Python...error.');
        end        
        %waitfor(errordlg(sprintf('%s',cmdout)));
        disp(cmdout);
        error('%s has not been installed properly.', ...
            upper(prgwkdir));
    end
    if isvalid(fw) 
        gui.gui_waitbar(fw,[],'Checking Python environment is complete.');
    end
    ok=true;

