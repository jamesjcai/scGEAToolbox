rng default
% a=rand(700);
% a=double(a>0.8);
% a=a'+a;
% a=a-diag(diag(a));
% g=graph(a);
% figure;
% plot(g)

addpath('C:\Users\jcai\Documents\GitHub\SBEToolbox_lite');
addpath('..\visualization_and_evaluation\');
%[x]=randnet_sw(10,1);
[x]=randnet_er(10,0.3);
[a1,p1,c1]=i_ceplot(x);

x(4,6)=false; x(6,4)=false;
x(4,2)=false; x(2,4)=false;
x(5,3)=true; x(3,5)=true;
[a2,~,c2]=i_ceplot(x,[],p1);
%[a3]=i_ceplot(x,a1(:,1));
% [coords_2_aligned, aligned_corr] = angular_alignment(a1(:,1), a2(:,1), 100);


function [coords_emb,p,coords]=i_ceplot(x,a1,p1)
    if nargin<3, p1=[]; end
    if nargin<2, a1=[]; end    
    coords_emb = coalescent_embedding(double(x), ...
        'RA1', 'LE', 'original', 2);
    if ~isempty(a1)
        [a2] = angular_alignment(a1, coords_emb(:,1), 100);
        coords_emb(:,1)=a2;
    end
    coords=coords_emb;
    [coords(:,1),coords(:,2)] = pol2cart(coords_emb(:,1),coords_emb(:,2));
    coords=coords./vecnorm(coords,2,2);
    figure; 
    subplot(1,2,1)
    g=graph(x);
    p=plot(g);
    %p.layout('auto');
    if ~isempty(p1)
        p.XData=p1.XData;
        p.YData=p1.YData;
    end
    subplot(1,2,2)
    scatter(coords(:,1),coords(:,2),[],g.centrality('degree'));
    text(coords(:,1),coords(:,2),string(1:size(coords,1)));
    t = linspace(0,2*pi,1000);
    line(sin(t),cos(t));
    axis equal
end

