function [done] = i_setnetwd(src, ~)

[parentfig] = gui.gui_getfigsce(src);
[done] = gui.i_setwrkdir('netanalywrkpath');
if done
     waitfor(gui.myHelpdlg(parentfig, ...
         "Network analysis working directory is " + ...
         "set successfully."));
end
