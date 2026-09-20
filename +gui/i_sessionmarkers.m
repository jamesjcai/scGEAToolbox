function [Tm, srcname] = i_sessionmarkers(newTm, newsrcname)
%I_SESSIONMARKERS Remember one customized marker gene list for the session.
%
%   [Tm, srcname] = gui.i_sessionmarkers()          read what is remembered
%   gui.i_sessionmarkers(Tm, srcname)               remember this list
%   gui.i_sessionmarkers([])                        forget it
%
%   Tm is the two-column table PKG.I_PARSEMARKERLIST returns - Var1 the cell
%   type name, Var2 its markers as one upper case comma separated string - and
%   SRCNAME a short description of where it came from ("typed list",
%   "markers.txt", ...) shown in the menus that offer to reuse it. Both are ""
%   / [] when nothing is remembered.
%
%   The list lives in a persistent variable, so it lasts as long as the MATLAB
%   session and is shared by every figure. Typing a marker list is slow enough
%   that retyping it for the next selection is the main reason the customized
%   path goes unused; nothing is written to disk, because a list that outlived
%   the session would then be a setting nobody remembers setting.
%
%   CLEAR FUNCTIONS forgets it, as does gui.i_sessionmarkers([]).
%
%   See also gui.i_getcustommarkers, pkg.i_parsemarkerlist,
%   gui.callback_Brush4Celltypes.

persistent storedTm storedsrcname

if nargin > 0
    if isempty(newTm)
        storedTm = [];
        storedsrcname = "";
    else
        storedTm = newTm;
        if nargin < 2 || strlength(string(newsrcname)) == 0
            storedsrcname = "customized list";
        else
            storedsrcname = string(newsrcname);
        end
    end
end

Tm = storedTm;
srcname = storedsrcname;
if isempty(srcname), srcname = ""; end
end
