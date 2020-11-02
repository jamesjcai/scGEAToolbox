function [] = plotADP(x,y,cut_indices,mass,I,cm_bool,rgb_c,linestyle)
% plots curve, data, to the current figure. If I and mass are specified as
% nonempty and cm_bool = 1 the curve of the projecting center of masses
% will also be plotted. Can leave F and iter as empty if movie is not being
% made.

[m,d] = size(y);
cut_partition = [0;cut_indices;m];


y_mass = zeros(m,1);
x_means = zeros(m,d);

if ~isempty(I) && cm_bool == 1
    assert(~isempty(mass) && length(mass) == length(I));
    for i=1:m
        y_mass(i) = sum(mass(I == i));
        x_means(i,:) = sum(bsxfun(@times,mass(I==i),x(I==i,:)'),2)/y_mass(i);
    end
end
singletons = [1;cut_indices+1];
singletons = singletons([1;cut_indices+1]==[cut_indices;m]);
if d==2
    if ~isempty(x)
        %scatter(x(:,1),x(:,2),15,[.7,.7,.7],'.');
        scatter(x(:,1),x(:,2),25,[0,0,0],'.');
        %scatter(x(:,1),x(:,2),20,[0.1,0.1,0.1],'filled');
%         if ~isempty(mass)
%             max_m = max(mass);
%             grey_scale = repmat((max_m-mass)'/max_m,1,3);
%             scatter(x(:,1),x(:,2),25,grey_scale,'.');
%         else
%             scatter(x(:,1),x(:,2),25,[0,0,0],'.');
%         end
    end
    
    hold on;
    for j=2:length(cut_indices)+2
        comp_y = y(cut_partition(j-1)+1:cut_partition(j),:);
        if length(comp_y(:,1)) > 1
            plot(comp_y(:,1),comp_y(:,2),'Color',rgb_c,'LineWidth',3,'LineStyle',linestyle);
            %scatter(comp_y(:,1),comp_y(:,2),4,'red');
        end
        if cm_bool == 1
            comp_x_means = x_means(cut_partition(j-1)+1:cut_partition(j),:);
            plot(comp_x_means(:,1),comp_x_means(:,2),'Color', 'cyan', 'LineWidth',2);
            scatter(comp_x_means(:,1),comp_x_means(:,2),4,'magenta');
        end
    end
    scatter(y(singletons,1),y(singletons,2),30,rgb_c,'filled','o');
    
else if d==3
        plot3dADP(x,y,cut_indices,cut_partition,x_means,cm_bool,rgb_c,linestyle);
        scatter3(y(singletons,1),y(singletons,2),y(singletons,3),30,rgb_c,'filled','o');
    else % if d>3, visualize in first 3 principal components
        n = length(x(:,1));
        if n-1>=3
            mean_x = mean(x,1);
            coeff = pca(x-repmat(mean_x,n,1));
            x_red = x*coeff(:,1:3);
            y_red = y*coeff(:,1:3);
            plot3dADP(x_red,y_red,cut_indices,cut_partition,x_means,cm_bool,rgb_c,linestyle);
            scatter3(y_red(singletons,1),y_red(singletons,2),y_red(singletons,3),30,rgb_c,'filled','o');
        end
    end
    grid off;
end
hold off;
axis equal;
end

function [] = plot3dADP(x,y,cut_indices,cut_partition,x_means,cm_bool,rgb_c,linestyle)

if ~isempty(x)
    scatter3(x(:,1),x(:,2),x(:,3),25,[0,0,0],'.');
end
hold on;
for j=2:length(cut_indices)+2
    comp_y = y(cut_partition(j-1)+1:cut_partition(j),:);
    if length(comp_y(:,1)) > 1
        plot3(comp_y(:,1),comp_y(:,2),comp_y(:,3),'Color',rgb_c,'LineWidth',2,'LineStyle',linestyle);
    end
    if cm_bool == 1
        comp_x_means = x_means(cut_partition(j-1)+1:cut_partition(j),:);
        plot3(comp_x_means(:,1),comp_x_means(:,2),comp_x_means(:,3),'Color', 'cyan', 'LineWidth',2);
        scatter3(comp_x_means(:,1),comp_x_means(:,2),comp_x_means(:,3),4,'magenta');
    end
end
end

