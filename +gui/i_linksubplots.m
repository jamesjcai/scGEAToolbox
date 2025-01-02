function i_linksubplots(~, ~)
%evalin('base', 'h=findobj(gcf,''type'',''axes'');');
%evalin('base', 'hlink = linkprop(h,{''CameraPosition'',''CameraUpVector''});');
%evalin('base', 'rotate3d on');

% evalin('base', 'linkprop(findobj(gcf,''type'',''axes''), {''CameraPosition'',''CameraUpVector''});');

% The linkprop function allows for more granular control by linking 
% specific properties across graphics objects. This is useful when you 
% want to synchronize properties other than axis limits or when linkaxes 
% doesn't provide the needed flexibility.
% % Link the x-axes of both subplots
% linkaxes([ax1, ax2], 'x');

h=gcf;
hlink = linkprop(findobj(h,'type','axes'), {'CameraPosition','CameraUpVector'});
setappdata(h,'UserData',hlink);
% When using linkprop, maintain a reference to the link object (hlink) in 
% the appropriate scope to prevent it from being cleared by MATLAB's 
% garbage collection. This can be done by storing it in the UserData 
% property of a figure or using setappdata. 

end