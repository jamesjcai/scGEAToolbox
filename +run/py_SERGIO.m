function [X]=py_SERGIO(A,ncells)

arguments    
    A (:,:) {mustBeNumericOrLogical,mustBeSquare(A)} = randnet_example
    ncells (1,1) {mustBePositive, mustBeInteger} = 1000
end
%ngenes (1,1) {mustBePositive, mustBeInteger} = 5

    ngenes=size(A,1);
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

    %tmpfilelist={'input.mat','output.mat','regs.txt','targets.txt'};
    tmpfilelist={'input.mat','output.mat'};
    if ~isdebug, pkg.i_deletefiles(tmpfilelist); end

    pkg_e_writesergiogrn(A);

   
    % fw=gui.gui_waitbar([],[],'Running GenKI...');
    save('input.mat','-v7.3','ncells','ngenes');
    
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


function mustBeSquare(X)
    if ~isequal(size(X,1),size(X,2)) || ~ismatrix(X)
        eid = 'Type:notSquareMatrix';
        msg = 'Must be a square matrix.';
        throwAsCaller(MException(eid,msg))
    end
end

function A=randnet_example
    rng(244)
    A=rand(5);
    A=A-diag(diag(A));
    A=A>0.55;
end

% function mustBeSquareEqualSize(X,n)
%     if ~isequal(size(X,1),n) || ~isequal(size(X,2),n)
%         eid = 'Size:notEqual';
%         msg = 'Size of A must equal number of genes.';
%         throwAsCaller(MException(eid,msg))
%     end
% end

% function mustBeRealUpperTriangular(a)
%     if ~(istriu(a) && isreal(a))
%         eidType = 'mustBeRealUpperTriangular:notRealUpperTriangular';
%         msgType = 'Input must be a real-valued, upper triangular matrix.';
%         throwAsCaller(MException(eidType,msgType))
%     end
% end
