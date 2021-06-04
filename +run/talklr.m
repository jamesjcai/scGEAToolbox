function [Tok,OUT]=talklr(X,g,c)

% https://doi.org/10.1101/2020.02.01.930602

if nargin<2, g=[]; c=[]; end
if isa(X, 'SingleCellExperiment')    
    g=X.g;
    c=X.c_cell_type_tx;
    X=X.X;
end
X=full(X);
[cx,cL]=grp2idx(c);
n=numel(cL);

pw=fileparts(mfilename('fullpath'));
dbfile=fullfile(pw,'..','resources','Ligand_Receptor.mat');
load(dbfile,'ligand','receptor','T');
T=T(:,2:6);

g=upper(g);

methodid=1;
switch methodid
    case 1
        X=sc_transform(X);
        %X=sc_norm(X);
        Xm=grpstats(X',cx,@mean)';
        cutoff=1;
    case 2
        Xm=zeros(numel(g),n);
        for k=1:n
            [t]=sc_hvg(X(:,cx==k),g,false);
            [y,id]=ismember(g,t.genes);
            Xm(y,k)=table2array(t(id(y),4));
        end
        cutoff=0;
end

idx=sum(Xm>cutoff,2)>0;
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
[n2]=size(ligand_mat,2);
assert(n==n2);

a1=ligand_mat(:,1).*receptor_mat;
a2=ligand_mat(:,2).*receptor_mat;
a3=ligand_mat(:,3).*receptor_mat;
M=[a1 a2 a3];
M=M./sum(M,2);
M=M.*log2(M*(n^2));
M(isnan(M))=0;
KL=real(sum(M,2));
Tok=[Tok, table(KL)];
[Tok,idx]=sortrows(Tok,'KL','descend');
%%
ligand_mat=ligand_mat(idx,:);
receptor_mat=receptor_mat(idx,:);
ligandok=ligandok(idx);
receptorok=receptorok(idx);
assert(isequal(size(ligand_mat,1),length(ligandok)))
KL=KL(idx);

OUT.cL=cL;
OUT.ligand_mat=ligand_mat;
OUT.receptor_mat=receptor_mat;
OUT.ligandok=ligandok;
OUT.receptorok=receptorok;
OUT.KL=KL;

for k=1:min([5 size(ligand_mat,1)])
    gui.i_crosstalkgraph(OUT,k);
end
end