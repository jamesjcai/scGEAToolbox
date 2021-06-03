function talkr(X,genelist,c_celltype_tx)

pw1=fileparts(mfilename('fullpath'));
pth=fullfile(pw1,'..','resources','Ligand_Receptor.mat');
load(pth)
%return;
%load Ligand_Receptor.mat

T=readtable('glom_normal_data.txt');
cutoff=4;
idx=sum(table2array(T(:,2:end))>cutoff,2)>0;
T=T(idx,:);

%%

ix=ismember(ligand,string(T.genes));
iy=ismember(receptor,string(T.genes));
idx=ix&iy;
ligandok=ligand(idx);
receptorok=receptor(idx);

[y1,idx1]=ismember(ligandok,string(T.genes));
[y2,idx2]=ismember(receptorok,string(T.genes));
assert(all(y1))
assert(all(y2))

t1=T(idx1,2:end);
t2=T(idx2,2:end);
t1.Properties.VariableNames={'a','b','c'};
t2.Properties.VariableNames={'d','e','f'};

Tx=[table(ligandok,receptorok),t1,t2];

%%
ligand_mat=table2array(t1)+1;
receptor_mat=table2array(t2)+1;
[p,n]=size(ligand_mat);
M=zeros(p,1);

a1=ligand_mat(:,1).*receptor_mat;
a2=ligand_mat(:,2).*receptor_mat;
a3=ligand_mat(:,3).*receptor_mat;
M=[a1 a2 a3];
M=M./sum(M,2);
M=M.*log2(M*(n^2));
M(isnan(M))=0

KL=sum(M,2);

%%
m=ligand_mat(1,:)'*receptor_mat(1,:);
m=m./sum(m(:));
g=digraph(m>0.3);
plot(g)

% final_mat=M(1,:)./sum(M(1,:))
