function [X]=e_uint2sparse(X)

if ~issparse(X) && ~isa(X,'double')
    try    
        X=sparse(double(X));
    catch ME
        if (strcmp(ME.identifier,'MATLAB:array:SizeLimitExceeded'))
            disp('Converting X to sparse.');
            % tic
            % S=spalloc(size(X,1),size(X,2),nnz(X));
            % idx=find(X>0);
            % S(idx)=X(idx);
            % toc
            % X=S;
            a=floor(size(X)./2);
            x1=sparse(double(X(1:a(1),1:a(2))));
            x2=sparse(double(X(a(1)+1:end,1:a(2))));
            x3=sparse(double(X(1:a(1),a(2)+1:end)));
            x4=sparse(double(X(a(1)+1:end,a(2)+1:end)));
            X=[x1 x3; x2 x4];
        end
    end
end
