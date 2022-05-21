function [T]=scInTime(X,genelist,ptime,varargin)

   p = inputParser;
   addOptional(p,'ncom',3,@(x) fix(x)==x & x>0);
   addOptional(p,'nsubsmpl',10,@(x) fix(x)==x & x>0);
   addOptional(p,'csubsmpl',500,@(x) fix(x)==x & x>0);
   addOptional(p,'savegrn',true,@islogical);
   addOptional(p,'mksparse',false,@islogical);
   parse(p,varargin{:});
   nsubsmpl=p.Results.nsubsmpl;
   csubsmpl=p.Results.csubsmpl;
   ncom=p.Results.ncom;
   savegrn=p.Results.savegrn;
   mksparse=p.Results.mksparse;

   [XM]=ten.i_nct(X,ptime,nsubsmpl,ncom,csubsmpl,savegrn,mksparse);
   assert(length(genelist)==size(XM,1))
   assert(size(X,1)==size(XM,1))
   n=length(genelist);
   R=reshape(corr(reshape(XM,[n.^2,10])',(1:10)'),[n,n]);
   R(isnan(R))=1;
   s=tsne(R,"NumDimensions",3);
   T.s=s;
   save('scInTime_res','T');




%%   
%    a=rand(3,3,10);
%    for k=1:10
%        a(:,:,k)=a(:,:,k)+k;
%    end
%    R=zeros(3,3);
%    for k=1:3
%        for l=1:3
%            R(k,l)=corr(reshape(a(k,l,:),10,1),[1:10]');
%        end
%    end
%    R2=reshape(corr(reshape(a,[9,10])',[1:10]'),3,3);
%    isequal(R,R2)




