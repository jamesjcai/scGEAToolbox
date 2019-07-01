function [X,Y]=run_combat(X,Y)
warning('run_combat is not recommended. Use run_combat2 instead.')
pw1=fileparts(which(mfilename));
addpath(fullfile(pw1,'thirdparty/ComBat'));
% addpath('thirdparty/ComBat/');
n1=size(X,2);
n2=size(Y,2);
batch=[ones(1,n1) 2*ones(1,n2)];
mod=ones(size(batch))';
normxy=combat([X Y],batch,mod);

% normxy(normxy<0)=0;
% normxy=round(normxy);

X=normxy(:,1:n1);
Y=normxy(:,n1+1:end);
