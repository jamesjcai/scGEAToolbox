function [sample_likelihoods,T]=py_MELD(X,batchid)
% MELD - a graph signal processing tool used to smooth a binary variable on 
% the cell-cell graph to determine which regions of its underlying 
% data manifold are enriched or depleted in cells with a specific 
% feature.
arguments    
    X (:,:) {mustBeNumeric}
    batchid (1,:) {mustBePositive, mustBeInteger}
end
sample_likelihoods=[];
T=[];

isdebug=false;
%if nargin<2
    %batchid=string([true(ceil(size(X,2)/2),1); false(floor(size(X,2)/2),1)]);
%    error('USAGE: score=run.MELD(X,batchid)');
%end
oldpth=pwd();
[pyok,wrkpth,x]=run.pycommon('py_MELD');
if ~pyok, return; end


%pw1=fileparts(mfilename('fullpath'));
%wrkpth=fullfile(pw1,'external','py_MELD');
%cd(wrkpth);

tmpfilelist={'batchid.txt','input.txt','output.txt','input.mat'};
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end

if issparse(X), X=full(X); end
X=sc_norm(X);
X=sqrt(X);

save('input.mat','X','batchid','-v7.3');
%save('input_v7.mat','X','batchid','-v7');
disp('Input file written.');
    
   
    fw=gui.gui_waitbar([],[],'Running ...');
    cmdlinestr=sprintf('"%s" "%s%sscript.py"', ...
        x.Executable,wrkpth,filesep);
    disp(cmdlinestr)
    [status]=system(cmdlinestr,'-echo');
    if isvalid(fw)
        gui.gui_waitbar(fw,[],'Running is complete');
    end
    if status==0 && exist('output.mat','file')
        load("output.mat","sout")        
    end

% x=pyenv;
% pkg.i_add_conda_python_path;
% if usematinput
%     if useh5
%         cmdlinestr=sprintf('"%s" "%s%sscript_h5.py"',x.Executable, ...
%             wrkpth,filesep);
%     else
%         cmdlinestr=sprintf('"%s" "%s%sscript_v7.py"',x.Executable, ...
%             wrkpth,filesep);
%     end
% else
%     cmdlinestr=sprintf('"%s" "%s%sscript_csv.py"',x.Executable, ...
%         wrkpth,filesep);
% end
% disp(cmdlinestr)
% [status]=system(cmdlinestr);
% 
% sample_likelihoods=[];
% if status==0 && exist('output.txt','file')
%     warning off
%     T=readtable('output.txt',"ReadVariableNames",true);
%     warning on
%     sample_likelihoods=table2array(T);
% end

if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
cd(oldpth);
end
