function [v]=GCL(X,k)

if nargin<2, k=50; end

pw1=fileparts(mfilename('fullpath'));
pth=fullfile(pw1,'thirdparty','GCL');
if ~(ismcc || isdeployed)
    addpath(pth);
end
%pth
[v] = gcl_ori(X, k);

end