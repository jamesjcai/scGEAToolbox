function [idx]=py_geosketch(X,n)

    isdebug = false;
    if nargin<2    
        n=min([2000,round(0.75*size(X,2))]);
    end
    idx=[];
    
    oldpth=pwd();
    pw1=fileparts(mfilename('fullpath'));
    wrkpth=fullfile(pw1,'external','py_geosketch');
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
        cd(oldpth);
        waitfor(errordlg(sprintf('%s',cmdout)));
        error('Python geosketch has not been installed properly.');
    end
    
    
    
    tmpfilelist={'input.mat','output.mat'};
    
    if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
    
    save('input.mat','-v7.3','X','n');
    disp('Input file written.');
    
    if isvalid(fw) 
        gui.gui_waitbar(fw,[],'Checking Python environment is complete');
    end
    
    fw=gui.gui_waitbar([],[],'Running geosketch...');
    cmdlinestr=sprintf('"%s" "%s%sscript.py"', ...
        x.Executable,wrkpth,filesep);
    disp(cmdlinestr)
    [status]=system(cmdlinestr,'-echo');
    if isvalid(fw)
        gui.gui_waitbar(fw,[],'Running geosketch is complete');
    end
    
    % if status==0 && exist('output.txt','file')
    %     T=readtable('output.txt');
    %     idx=T.Var1+1;
    % end

    if status==0 && exist('output.mat','file')
        load("output.mat","idx")
        idx=idx+1;
    end    
    
    if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
    cd(oldpth);
end
