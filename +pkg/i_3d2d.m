function [s2]=i_3d2d(s3,ax,bx)

A = viewmtx(ax,bx);
x4d = [s3,ones(size(s3,1),1)]';
x2d = A*x4d;
s2=x2d';
s2=s2(:,1:3);