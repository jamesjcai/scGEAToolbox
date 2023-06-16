function pkg_e_writesergiogrn(A)

if nargin<1
    rng(244)
    A=rand(5);
    A=A-diag(diag(A));
    A=A>0.55;
end

g=digraph(A);
% plot(g)
d=indegree(g);
A=full(adjacency(g));

% a1=triu(A);
% a2=tril(A);

fid1=fopen('targets.txt','w');
fid2=fopen('regs.txt','w');
for k=1:size(g.Nodes,1)
    if d(k)==0
        fprintf(fid2,'%d, 0.5\n',k-1);
    else
        %[d(k) sum(A(:,k))]
        f1=find(A(:,k));
        f2=2*ones(size(f1));
        f3=0.5*ones(size(f1));
        rd=ones(size(f1));
        rd(rand(size(f1))>0.8)=-1;
        f3=f3.*rd;
        s=sprintf('%d, %d, %s%s%s',k-1,d(k), ...
            sprintf('%d, ',f1-1), ...
            sprintf('%f, ',f3), ...
            sprintf('%d, ',f2));        
        s=strtrim(s);
        s=extractBefore(s,strlength(s));
        %disp(s);
        fprintf(fid1,'%s\n', s);
    end
end
fclose(fid1);
fclose(fid2);

