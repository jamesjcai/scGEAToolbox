function [X]=py_SERGIO   

    X=[];
    isdebug = false;
    oldpth=pwd();
    pw1=fileparts(mfilename('fullpath'));
    wrkpth=fullfile(pw1,'external','py_SERGIO');
    cd(wrkpth);   
    
    fw = gui.gui_waitbar([],[],'Checking Python environment...');    
    x=pyenv;
    try
        pkg.i_add_conda_python_path;
    catch
        
    end
    cmdlinestr=sprintf('"%s" "%s%srequire.py"', ...
            x.Executable,wrkpth,filesep);
    disp(cmdlinestr)
    [status,cmdout]=system(cmdlinestr,'-echo');
    if status~=0
        gui.gui_waitbar(fw,true);
        cd(oldpth);
        waitfor(errordlg(sprintf('%s',cmdout)));
        error('Python SERGIO has not been installed properly.');
    end 
    
    if isvalid(fw) 
        gui.gui_waitbar(fw,[],'Checking Python environment is complete');
    end

    tmpfilelist={'output.mat'};
    if ~isdebug, pkg.i_deletefiles(tmpfilelist); end

   
    % fw=gui.gui_waitbar([],[],'Running GenKI...');
    
    cmdlinestr=sprintf('"%s" "%s%sscript.py"', ...
        x.Executable,wrkpth,filesep);
    disp(cmdlinestr)
    [status]=system(cmdlinestr,'-echo');

    % if isvalid(fw)
    %     gui.gui_waitbar(fw,[],'Running GenKI is complete');
    % end
    
     if status==0 && exist('output.mat','file')
         load('output.mat','X');
     end

    if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
    cd(oldpth);
end
