function [s]=r_node2vec(A,g)

narginchk(1,2);
if isa(A,'graph') || isa(A,'digraph') 
    A=G.adjacency;
    g=G.Nodes.Name;
else
    if nargin<2
        g=string((1:size(A,1))');
    end
end

s=[];
isdebug=false;
oldpth=pwd();
[isok,msg]=commoncheck_R('R_node2vec');
if ~isok, error(msg); end

tmpfilelist={'input.txt','output.txt'};
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end

fid=fopen('input.txt','w');
fprintf(fid,'gene1,gene2\n');
n=size(A,1);
for i=1:n-1
for j=i+1:n
    if A(i,j)>0
        fprintf(fid,'%s,%s\n',g(i),g(j));
    end
end
end
fclose(fid);


Rpath=getpref('scgeatoolbox','rexecutablepath');
pkg.RunRcode('script.R',Rpath);

if ~exist('output.txt','file'), return; end
T=readtable('output.txt');
g2=string(T.Var1);
s=table2array(T(:,2:end));
[y,idx]=ismember(g,g2);
assert(all(y));
assert(isequal(g(:),g2(idx)));
s=s(idx,:);
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
cd(oldpth);
end
