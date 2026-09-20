function [capfig, captab, captitle] = i_pptxslideplan(F, glist, tabs)
%I_PPTXSLIDEPLAN Work out what each slide of a PowerPoint export shows.
%
%   [capfig, captab, captitle] = gui.i_pptxslideplan(F, glist, tabs)
%
% Inputs:
%   F      cell array of figure handles
%   glist  slide titles to go with them, {''} or {} for untitled slides
%   tabs   tabs of F{1} to expand into a slide each, empty for none
%
% Outputs, one element per slide:
%   capfig   the figure to photograph
%   captab   the tab to bring forward first, [] when there is none
%   captitle the slide title, "" for an untitled slide
%
% Two shapes, and the caller decides which by what it passes as TABS. With no
% tabs it is the original behaviour: a slide per figure, titled from GLIST.
% With tabs it is a slide per tab of the first figure, titled from the tab -
% which is how a gui.MYFIGURE carrying a tab per gene reaches the deck as
% more than the one gene that happened to be on top.
%
% Pure on purpose. The question of whether the user wants every tab belongs
% to GUI.I_EXPORT2PPTX, which has a figure to put a dialog on; keeping it out
% of here is what lets the slide list be checked without one.
%
% see also: gui.i_export2pptx, gui.i_figtabs

if nargin < 3, tabs = gobjects(0); end
if nargin < 2, glist = {}; end

capfig = {};
captab = {};
captitle = strings(0, 1);
if isempty(F), return; end

if ~isempty(tabs)
    tabs = tabs(:);
    capfig = repmat(F(1), numel(tabs), 1);
    captab = num2cell(tabs);
    % GET over a handle vector returns a cell column, one Title each; over a
    % single handle it returns the char itself, which STRING would otherwise
    % split into one title per character.
    captitle = reshape(string(get(tabs, 'Title')), [], 1);
    return;
end

capfig = F(:);
captab = repmat({[]}, numel(capfig), 1);
captitle = strings(numel(capfig), 1);

% {''} is what GUI.I_SAVEMAINFIG passes for "no titles", so an empty first
% entry means the whole list is decoration rather than a title per slide.
if ~isempty(glist) && ~isempty(glist{1})
    n = min(numel(glist), numel(capfig));
    captitle(1:n) = reshape(string(glist(1:n)), [], 1);
end
end
