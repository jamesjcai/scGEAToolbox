function [OUT]=talklr(X,g,c)
%TALKLR - uncovers ligand-receptor mediated intercellular crosstalk
%
% Usage: [OUT]=run.talklr(X,g,c)
% X - gene-by-cell expression matrix
% g - list of genes (string array)
% c - index of cell types
%
% This function implements the method of Yuliang Wang (2020)
% https://doi.org/10.1101/2020.02.01.930602

if nargin<2, g=[]; c=[]; end
if isa(X, 'SingleCellExperiment')    
    g=X.g;
    c=X.c_cell_type_tx;
    X=X.X;
end


[M,OUT]=ii_talkr(X,g,c);

[n]=size(OUT.ligand_mat,2);
M=M.*log2(M*(n^2));
M(isnan(M))=0;
KL=real(sum(M,2));
Tok=[OUT.Tok, table(KL)];
[Tok,idx]=sortrows(Tok,'KL','descend');
%%
OUT.ligand_mat=OUT.ligand_mat(idx,:);
OUT.receptor_mat=OUT.receptor_mat(idx,:);
OUT.ligandok=OUT.ligandok(idx);
OUT.receptorok=OUT.receptorok(idx);
% assert(isequal(size(ligand_mat,1),length(ligandok)))
OUT.KL=KL(idx);
OUT.T=Tok;

disp('Try: gui.i_crosstalkgraph(OUT,k);')
end


function [M,OUT]=ii_talkr(X,g,c)
    if isa(X,'SingleCellExperiment')    
        c=X.c_cell_type_tx;
        g=X.g;
        X=X.X;
    end

    [cx,cL]=grp2idx(c);
    n=numel(cL);

    pw=fileparts(mfilename('fullpath'));
    dbfile=fullfile(pw,'..','resources','Ligand_Receptor.mat');
    load(dbfile,'ligand','receptor','T');
    %T=T(:,2:6);

    g=upper(g);

    methodid=1;
    switch methodid
        case 1
            X=sc_transform(X,'type','kNNSmoothing');
            % X=sc_transform(X)
            % X=sc_norm(X);
            % Xm=grpstats(X',cx,@mean)';
            [Xm]=i_grpmean(X,cx);
            cutoff=1.0;
        case 2
            Xm=zeros(numel(g),n);
            for k=1:n
                [t]=sc_hvg(X(:,cx==k),g,false);
                [y,id]=ismember(g,t.genes);
                Xm(y,k)=table2array(t(id(y),4));
            end
            Xm(Xm<0)=0;
            cutoff=0.05;
    end

    idx=sum(Xm>cutoff,2)>1;
    Xm=Xm(idx,:);
    g=g(idx);

    %%

    ix=ismember(ligand,g);
    iy=ismember(receptor,g);
    idx=ix&iy;
    ligandok=ligand(idx);
    receptorok=receptor(idx);
    Tok=T(idx,:);

    [y1,idx1]=ismember(ligandok,g);
    [y2,idx2]=ismember(receptorok,g);
    assert(all(y1))
    assert(all(y2))

    %ligandok=ligandok(idx1);
    %receptorok=receptorok(idx2);

    ligand_mat=Xm(idx1,:);
    receptor_mat=Xm(idx2,:);
    % t1.Properties.VariableNames={'a','b','c'};
    % t2.Properties.VariableNames={'d','e','f'};
    % Tx=[table(ligandok,receptorok),t1,t2];

    %%
    [p,n2]=size(ligand_mat);
    assert(n==n2);


    M=zeros(p,n.^2);
    for k=1:n
        M(:,(n*(k-1)+1):n*k)=ligand_mat(:,k).*receptor_mat;
    end
    M=M./sum(M,2);
    %M=M.*log2(M*(n^2));
    %M(isnan(M))=0;

    %M=zeros(size(ligand_mat,1),n^2);
    %M=[];
    %for k=1:n
    %    M=[M ligand_mat(:,k).*receptor_mat];
    %end
    % a2=ligand_mat(:,2).*receptor_mat;
    % a3=ligand_mat(:,3).*receptor_mat;
    % M=[a1 a2 a3];
    %M=M./sum(M,2);
    %size(M)

    OUT.cL=cL;
    OUT.ligand_mat=ligand_mat;
    OUT.receptor_mat=receptor_mat;
    OUT.ligandok=ligandok;
    OUT.receptorok=receptorok;
    OUT.Tok=Tok;
end


