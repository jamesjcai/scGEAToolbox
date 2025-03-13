function [done] = i_setextwd(src, ~)
%I_SETEXTWD - set externel workding directory

%see also: I_SETPYENV, I_SETRENV 
[parentfig, ~] = gui.gui_getfigsce(src);
[done] = gui.i_setwrkdir('externalwrkpath');
if done
     gui.myHelpdlg(parentfig, "External program working " + ...
         "directory is set successfully.");
end

