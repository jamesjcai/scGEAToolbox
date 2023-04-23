function [X,genelist]=sc_qcfilter_lite(X,genelist,libszcutoff,...
                    mtratio,min_cells_nonzero,gnnumcutoff)
if nargin<6, gnnumcutoff=200; end
if nargin<5 || isempty(min_cells_nonzero), min_cells_nonzero=10; end
if nargin<4 || isempty(mtratio), mtratio=0.15; end          % 0.10
if nargin<3 || isempty(libszcutoff), libszcutoff=500; end   % 1000


idx=startsWith(genelist,"mt-",'IgnoreCase',true);
if sum(idx)>0
    lbsz=sum(X,1);
    lbsz_mt=sum(X(idx,:),1);
    f_mtreads=lbsz_mt./lbsz;
    keptidx=f_mtreads<mtratio;
    if sum(~keptidx)>0
        X(:,~keptidx)=[];
    end
    
end


oldsz=0;
newsz=1;

while ~isequal(oldsz,newsz)
    oldsz=size(X);
    % [X,genelist]=sc_filterg(X,genelist);   % remove empty genes
    [u]=sum(X,2);
    idx=u<1;
    X(idx,:)=[];
    genelist=genelist(~idx);
    

    % [X,keptidx]=sc_filterc(X);             % remove empty cells
    lbsz=sum(X,1);
    idx=lbsz<1;
    if ~isempty(idx), X(:,idx)=[]; end


    %[X,genelist]=sc_selectg(X,genelist,min_cells_nonzero);
    nc=sum(X>=1,2);
    if min_cells_nonzero<1
        idx=nc>=min_cells_nonzero*size(X,2);
    else
        idx=nc>=min_cells_nonzero;
    end
    X(~idx,:)=[];
    genelist=genelist(idx);
    

    %[X,genelist]=sc_rmdugenes(X,genelist);
    %[X,keptidx]=sc_selectc(X,libszcutoff,gnnumcutoff);
    
    libsz=sum(X,1);
    gnnum=sum(X>0,1);
    
    if libszcutoff>1.0
        keptidx = (libsz>=libszcutoff) & (gnnum>=gnnumcutoff);
    else
        keptidx = (libsz>=quantile(libsz,libszcutoff)) & (gnnum>=gnnumcutoff);
    end
    X(:,~keptidx)=[];
    
    newsz=size(X);
end
end
