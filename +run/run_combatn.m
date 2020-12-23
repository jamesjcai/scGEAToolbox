function [X]=run_combatn(X,batchid)
pw1=fileparts(mfilename('fullpath'));
addpath(fullfile(pw1,'thirdparty/ComBat'));
if size(batchid,1)~=1
    batchid=batchid';
end
mod=ones(size(batchid))';
X=combat(X,batchid,mod);

