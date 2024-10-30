function i_linksubplots(~, ~)
%evalin('base', 'h=findobj(gcf,''type'',''axes'');');
%evalin('base', 'hlink = linkprop(h,{''CameraPosition'',''CameraUpVector''});');
%evalin('base', 'rotate3d on');
evalin('base', 'linkprop(findobj(gcf,''type'',''axes''), {''CameraPosition'',''CameraUpVector''});');
end