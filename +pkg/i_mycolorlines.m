function a=i_mycolorlines(kc)

% see also: pkg.distinguishable_colors

if kc <= 7 && kc>0
    a=lines(kc);
elseif kc>7 && kc<=12
%    colormap(gui.linspecer(kc,'qualitative'));    
%    colormap default;
     a=turbo(kc);
else
    a=turbo(kc);
    % colormap(gui.linspecer(kc,'sequential'));
    % colormap(gui.distinguishable_colors(kc));
    % colormap(pkg.i_mycolormap(kc));
    %
    % see also: i_gscatter3
    % see also: sc_scatter_sce;
    % see also: gui.sc_multigroupings
end

% if kc<=7 && kc>0
%     colormap(lines(kc));    
% else
% %     cx=colormap('autumn');
% %     cx(1,:)=[.8 .8 .8];
% %     colormap(cx);
%     % colormap default
%     colormap(gui.linspecer(kc));
% % https://www.mathworks.com/help/matlab/ref/colormap.html
% end


% a = pkg.distinguishable_colors(kc);
