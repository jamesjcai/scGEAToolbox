G=graph;
for k=1:size(Graph,2)
    g=Graph(:,k);
    for i=1:length(g)
        if g(i)~=k
            G=addedge(G,k,g(i));
        end
    end
end

%%
figure;
% G=simplify(G);
p=plot(G);

% p.XData=s(:,1)';
% p.YData=s(:,2)';
% p.ZData=s(:,3)';

[T,pred] = minspantree(G);
highlight(p,T,'EdgeColor','r')


%layout(p,'force3')
%view(3)







