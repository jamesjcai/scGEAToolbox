function [X,M]=run_mcimpute(X,donorm)

if nargin<2, donorm=true; end

pw1=fileparts(mfilename('fullpath'));
pth=fullfile(pw1,'thirdparty/McImpute');
addpath(pth);

if ~isnormalized(X) && ~donorm
    warning('Normalized X as input is recommended.');
end

libsize=sum(X)';

if donorm
    X=sc_norm(X,'type','libsize');
end

X=log2(X+1);
% McImpute needs [cells x genes]
X=X';

IDX = find(X);
M = opRestriction(numel(X),IDX);
y = M(X(:),1);
 
[Xrec] = IST_MC(y,M,size(X),0); %mask changed and lansvd changed, with NN constraint
%Xrec = IST_eMC(y,M,size(X),11); 

normed_data=2.^Xrec-1;
raw_data=(normed_data.*libsize)./ median(libsize);
M=round(raw_data);

X=normed_data';

