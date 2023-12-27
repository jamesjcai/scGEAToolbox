function [done] = i_setnetwd(~, ~)
[done] = gui.i_setwrkdir('netanalywrkpath');
if done
     waitfor(helpdlg("Network analysis working directory is set successfully.", ''));
end
