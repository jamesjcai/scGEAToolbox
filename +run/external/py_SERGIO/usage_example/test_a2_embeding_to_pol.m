D=squareform(pdist([s0; s1]));
[Y] = Isomap(D, 'k', 3);
s=Y.coords{2}';
coords_emb=s;
[coords_emb(:,1),coords_emb(:,2)] = cart2pol(s(:,1),s(:,2));


[coords(:,1),coords(:,2)] = pol2cart(coords_emb(:,1),coords_emb(:,2));
 coords=coords./vecnorm(coords,2,2);

    figure; 
    subplot(1,2,1)
        scatter(s(:,1),s(:,2))
        text(s(:,1),s(:,2),string(1:20))
        hold on
        for k=1:10
            if ismember(k,[2 4 6 3 5])
                line([s(k,1) s(k+10,1)],[s(k,2) s(k+10,2)],'color','r');
            else    
                line([s(k,1) s(k+10,1)],[s(k,2) s(k+10,2)])
            end
        end
    subplot(1,2,2)
        scatter(coords(:,1),coords(:,2));
        text(coords(:,1),coords(:,2),string(1:size(coords,1)));
        t = linspace(0,2*pi,1000);
        line(sin(t),cos(t));
        axis equal


        a1=coords_emb(1:10,:);
        a2=coords_emb(11:end,:);
        