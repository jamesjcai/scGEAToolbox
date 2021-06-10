function ndm=CSN_transform(X)
    pw1=fileparts(mfilename('fullpath'));
    pth=fullfile(pw1,'thirdparty','CSN_transform');
    addpath(pth);
    % G=csnet(X(:,1:10));
    % csnedge(X(1,:),X(2,:));
    ndm=csndm(X);
end