function s = i_classlistdivider()
%I_CLASSLISTDIVIDER  The row separating built-in groupings from attributes.
%
%   s = i_classlistdivider()
%
% A uilistbox has no per-item disable, so the separator has to be an item
% like any other and the pickers have to refuse it. Both GUI.I_SELECT1CLASS
% and GUI.I_SELECTNCLASS drop it from a selection, which makes choosing it
% alone behave as cancelling does; this function is the one place the text
% is written so those checks cannot drift from what is displayed.
%
% Built from CHAR rather than a literal: the glyph is U+2500 BOX DRAWINGS
% LIGHT HORIZONTAL, and building it by code point keeps every .m file in
% this package pure ASCII, so nothing depends on how a reader decodes them.
%
% See also GUI.I_SELECT1CLASS, GUI.I_SELECTNCLASS

rule = repmat(char(9472), 1, 6);
s = [rule ' Cell Attributes ' rule];
end
