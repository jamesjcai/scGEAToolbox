function [A]=sc_pcnetdenoised(X,varargin)
%Construct network using scTenifoldNet (accurate, robust, but slow)
% A=sc_pcnetdenoised(X);
%
% X is gene x cell matrix
% 

% if exist('sctenifoldnet','file')~=2
%     error('Requires sctenifoldnet.m');
% end

import ten.*

    if nargin<1
        error(sprintf('USAGE: A=ten.sc_pcnetdenoised(X);\n       A=sc_pcnetdenoised(X,''smplmethod'',''Jackknife'');'));
    end
    
   p = inputParser;
   addOptional(p,'smplmethod',"Jackknife",@(x) (isstring(x)|ischar(x))&ismember(lower(string(x)),["jackknife","bootstrap"]));
   addOptional(p,'tdmethod',"CP",@(x) (isstring(x)|ischar(x))&ismember(upper(string(x)),["CP","TUCKER"]));
   addOptional(p,'nsubsmpl',10,@(x) fix(x)==x & x>0);
   addOptional(p,'csubsmpl',500,@(x) fix(x)==x & x>0);
   addOptional(p,'savegrn',false,@islogical);
   addOptional(p,'donorm',true,@islogical);
   parse(p,varargin{:});
   tdmethod=p.Results.tdmethod;
   nsubsmpl=p.Results.nsubsmpl;
   csubsmpl=p.Results.csubsmpl;
   smplmethod=p.Results.smplmethod;
   savegrn=p.Results.savegrn;
   donorm=p.Results.donorm;
   
   switch upper(tdmethod)
       case "CP"
           tdmethod=1;
       case "TUCKER"
           tdmethod=2;
   end
   switch lower(smplmethod)
       case "jackknife"
           usebootstrp=false;
       case "bootstrap"
           usebootstrp=true;
   end

    
    if exist('@tensor/tensor.m','file')~=2
        error('Need Tensor Toolbox for MATLAB (https://www.tensortoolbox.org/)');
    end
    if exist('sc_pcnetpar.m','file')~=2
        error('Need sc_pcnetpar.m in scGEAToolbox https://github.com/jamesjcai/scGEAToolbox');
    end
    
    if donorm
        X=sc_norm(X,"type","libsize");
        X=log(1+X);
        %X=sc_transform(X);
    end
    
    %tic
    [XM]=ten.i_nc(X,nsubsmpl,3,csubsmpl,usebootstrp);
    %toc

    %tic
    disp('Tensor decomposition')
    [A]=ten.i_td1(XM,tdmethod);
    %toc
    if savegrn
        tstr=matlab.lang.makeValidName(string(datetime));
        save(sprintf('A_%s',tstr),'A','-v7.3');
    end    
end