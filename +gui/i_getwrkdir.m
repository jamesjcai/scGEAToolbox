function [wkdir] = i_getwrkdir(preftagname)
%I_GETWRKDIR - get workding directory
%see also: I_SETPYENV, I_SETRENV 
wkdir = [];
% preftagname = 'externalwrkpath';    i_setextwd
% preftagname = 'netanalywrkpath';    i_setnetwd
if nargin < 1, preftagname = 'externalwrkpath'; end
if ~ispref('scgeatoolbox', preftagname), return; end
wkdir = getpref('scgeatoolbox', preftagname, []);
if isempty(wkdir)
    [done] = gui.i_setwrkdir;
    if done
        wkdir = getpref('scgeatoolbox', preftagname, []);        
    end
end
