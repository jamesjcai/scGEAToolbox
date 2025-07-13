%% Plot the quantum bit state (either mixed state or pure state) in a
%%% Bloch sphere
%%% Written by Guoqing Wang gq_wang@mit.edu in Feb-2020, modified by Guoqing in
%%% Mar-2021
%%% Input "state" should be unit cell

%{
qtm.BlochSphere_Evolution_GW21({[0,1]',[1,1]'/sqrt(2), [1,0]'});
%}

function [figure1,x,y,z] = BlochSphere_Evolution_GW21(state,fig)
if nargin<2
    figure1 = figure;
    [x,y,z] = sphere();
    hold on;box on;set(gcf,'position',[700 100 800 600]);
    mesh(x,y,z,'FaceLighting','none','EdgeLighting','flat','FaceAlpha',0.1,'FaceColor',[0.83 0.81 0.78],'EdgeColor',[0.49 0.49 0.49]);
    alpha(0.1);
else
    figure1=fig;
end
curve_nbs = length(state);
x={};y={};z={};
for jj = 1:curve_nbs
    [m,n] = size(state{jj});
    if m==2 %% unitary, need to convert to xyz axis
        [x{jj},y{jj},z{jj}] = Blochsphere_Convertion_GW21(state{jj});
    else %%m=3, already xyz axis
        x{jj} = state{jj}(1,:);
        y{jj} = state{jj}(2,:);
        z{jj} = state{jj}(3,:);
    end
    plot3(x{jj},y{jj},z{jj},'o','linewidth',1.8);  % Plot
end
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
zlabel('z','FontSize',14);
% Example of using quiver3 correctly
%quiver3(0, 0, 0, 1, 0, 0, 1.2, 'LineWidth', 2, 'Color', 'k', 'Marker', 'o');

quiver3(0,0,0,1,0,0,1.2,'LineWidth',2,'color','k','Marker', 'o');
quiver3(0,0,0,0,1,0,1.2,'LineWidth',2,'color','k','Marker', 'o');
quiver3(0,0,0,0,0,1,1.2,'LineWidth',2,'color','k','Marker', 'o');
text(1.25,0,0,'x','FontSize',18);
text(0,1.25,0,'y','FontSize',18);
text(0,0,1.25,'z','FontSize',18);
title('Bloch sphere visualization');
view([116.74 16.56]);
axis equal;
end
%% This function convert a 2-by-1 qubit state to its 3D xyz coordinates
function [x,y,z,aa_phi,aa_theta] = Blochsphere_Convertion_GW21(aa)
    aa_phi = angle(aa(2,:))-angle(aa(1,:));
    aa_theta = 2*acos(abs(aa(1,:)));
    x = sin(aa_theta).*cos(aa_phi);
    y = sin(aa_theta).*sin(aa_phi);
    z = cos(aa_theta);
end