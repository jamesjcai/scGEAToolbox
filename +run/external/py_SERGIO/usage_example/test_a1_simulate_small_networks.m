rng default
addpath('C:\Users\jcai\Documents\GitHub\SBEToolbox_lite');
addpath('..\visualization_and_evaluation\');
%[x]=randnet_sw(10,1);
[A0]=randnet_er(10,0.3);
x=A0;
x(4,6)=false; x(6,4)=false;
x(4,2)=false; x(2,4)=false;
x(5,3)=true; x(3,5)=true;
A1=x;

    figure; 
    subplot(1,2,1)
        g0=graph(A0);
        p0=plot(g0);
    subplot(1,2,2)
        g1=graph(A1);
        p1=plot(g1);

g=string(1:10)';
s0=run.r_node2vec(A0,g);
s1=run.r_node2vec(A1,g);

%%
[~,s]=pca([s0;s1],'NumComponents',2);
%s=tsne([s0;s1],'NumDimensions',2);
figure;
scatter(s(:,1),s(:,2))
text(s(:,1),s(:,2),string(1:20))
%text(s(11:end,1),s(11:end,2),string(100*(1:10)))
hold on
for k=1:10
    if ismember(k,[2 4 6 3 5])
        line([s(k,1) s(k+10,1)],[s(k,2) s(k+10,2)],'color','r');
    else    
        line([s(k,1) s(k+10,1)],[s(k,2) s(k+10,2)])
    end
end


    


