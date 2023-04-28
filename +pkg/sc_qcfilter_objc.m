function [genelist,b]=sc_qcfilter_objc(Xobj,genelist,libszcutoff,...
                    mtratio,min_cells_nonzero,gnnumcutoff)
if nargin<6, gnnumcutoff=200; end
if nargin<5 || isempty(min_cells_nonzero), min_cells_nonzero=10; end
if nargin<4 || isempty(mtratio), mtratio=0.15; end          % 0.10
if nargin<3 || isempty(libszcutoff), libszcutoff=500; end   % 1000


idx=startsWith(genelist,"mt-",'IgnoreCase',true);
if sum(idx)>0
    lbsz=sum(Xobj.data,1);
    lbsz_mt=sum(Xobj.data(idx,:),1);
    f_mtreads=lbsz_mt./lbsz;
    keptidx=f_mtreads<mtratio;
    if sum(~keptidx)>0
        Xobj.data(:,~keptidx)=[];
    end
    
end


oldsz=0;
newsz=1;

while ~isequal(oldsz,newsz)
    oldsz=size(Xobj.data);
    % [X,genelist]=sc_filterg(X,genelist);   % remove empty genes
    [u]=sum(Xobj.data,2);
    idx=u<1;
    Xobj.data(idx,:)=[];
    genelist=genelist(~idx);
    

    % [X,keptidx]=sc_filterc(X);             % remove empty cells
    lbsz=sum(Xobj.data,1);
    idx=lbsz<1;
    if ~isempty(idx), Xobj.data(:,idx)=[]; end


    %[X,genelist]=sc_selectg(X,genelist,min_cells_nonzero);
    nc=sum(Xobj.data>=1,2);
    if min_cells_nonzero<1
        idx=nc>=min_cells_nonzero*size(Xobj.data,2);
    else
        idx=nc>=min_cells_nonzero;
    end
    Xobj.data(~idx,:)=[];
    genelist=genelist(idx);
    

    %[X,genelist]=sc_rmdugenes(X,genelist);
    %[X,keptidx]=sc_selectc(X,libszcutoff,gnnumcutoff);
    
    libsz=sum(Xobj.data,1);
    gnnum=sum(Xobj.data>0,1);
    
    if libszcutoff>1.0
        keptidx = (libsz>=libszcutoff) & (gnnum>=gnnumcutoff);
    else
        keptidx = (libsz>=quantile(libsz,libszcutoff)) & (gnnum>=gnnumcutoff);
    end
    Xobj.data(:,~keptidx)=[];
usr = memory;
b = usr.MemUsedMATLAB/1e6;

    newsz=size(Xobj.data);
end
end
