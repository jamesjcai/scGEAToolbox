function [X,Y]=mt_ComBat(X,Y)
pw1=fileparts(mfilename('fullpath'));
if ~(ismcc || isdeployed)
    addpath(fullfile(pw1,'external','mt_ComBat'));
end

n1=size(X,2);
n2=size(Y,2);
batch=[ones(1,n1) 2*ones(1,n2)];
mod=ones(size(batch))';
normxy=combat([X Y],batch,mod);

% normxy(normxy<0)=0;
% normxy=round(normxy);

X=normxy(:,1:n1);
Y=normxy(:,n1+1:end);
end