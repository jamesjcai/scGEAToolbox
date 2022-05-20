function [XM]=scInTime(X,genelist,ptime,varargin)

   p = inputParser;
   addOptional(p,'ncom',3,@(x) fix(x)==x & x>0);
   addOptional(p,'nsubsmpl',10,@(x) fix(x)==x & x>0);
   addOptional(p,'csubsmpl',500,@(x) fix(x)==x & x>0);
   addOptional(p,'savegrn',false,@islogical);
   parse(p,varargin{:});
   nsubsmpl=p.Results.nsubsmpl;
   csubsmpl=p.Results.csubsmpl;
   ncom=p.Results.ncom;
   savegrn=p.Results.savegrn;

   [XM]=ten.i_nct(X,ptime,nsubsmpl,ncom,csubsmpl,savegrn);
   a=rand(3,3,10);
   for k=1:10
       a(:,:,k)=a(:,:,k)+k;
   end
   R=zeros(3,3);
   for k=1:3
       for l=1:3
           R(k,l)=corr(reshape(a(k,l,:),10,1),[1:10]');
       end
   end
   R2=reshape(corr(reshape(a,10,9),[1:10]'),3,3);

   corr(reshape(a,10,9),[1:10]')

