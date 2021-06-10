function [M1,M2,OUT1,OUT2]=i_talkr2(X,g,c,b)
if isa(X,'SingleCellExperiment')
    b=X.c_batch_id;
    c=X.c_cell_type_tx;
    g=X.g;
    X=X.X;
end
[b]=grp2idx(b);
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
        Xm1=i_grpmean(X(:,b==1),cx(b==1));
        Xm2=i_grpmean(X(:,b==2),cx(b==2));
        
        %Xm1=grpstats(X(:,b==1)',cx(b==1),@mean)';
        %Xm2=grpstats(X(:,b==2)',cx(b==2),@mean)';
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

idx1=sum(Xm1>cutoff,2)>0;
idx2=sum(Xm2>cutoff,2)>0;
idx=idx1&idx2;

Xm1=Xm1(idx,:);
Xm2=Xm2(idx,:);
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

ligand_mat1=Xm1(idx1,:);
receptor_mat1=Xm1(idx2,:);
ligand_mat2=Xm2(idx1,:);
receptor_mat2=Xm2(idx2,:);
% t1.Properties.VariableNames={'a','b','c'};
% t2.Properties.VariableNames={'d','e','f'};
% Tx=[table(ligandok,receptorok),t1,t2];

%%
[n1]=size(ligand_mat1,2);
[n2]=size(ligand_mat2,2);
assert(n1==n2);


a1=ligand_mat1(:,1).*receptor_mat1;
a2=ligand_mat1(:,2).*receptor_mat1;
a3=ligand_mat1(:,3).*receptor_mat1;
M1=[a1 a2 a3];
M1=M1./sum(M1,2);


a1=ligand_mat2(:,1).*receptor_mat2;
a2=ligand_mat2(:,2).*receptor_mat2;
a3=ligand_mat2(:,3).*receptor_mat2;
M2=[a1 a2 a3];
M2=M2./sum(M2,2);


OUT1.cL=cL;
OUT1.ligand_mat=ligand_mat1;
OUT1.receptor_mat=receptor_mat1;
OUT1.ligandok=ligandok;
OUT1.receptorok=receptorok;
OUT1.Tok=Tok;

OUT2.cL=cL;
OUT2.ligand_mat=ligand_mat2;
OUT2.receptor_mat=receptor_mat2;
OUT2.ligandok=ligandok;
OUT2.receptorok=receptorok;
OUT2.Tok=Tok;

end


