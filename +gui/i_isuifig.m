function y = i_isuifig(parentfig)
% https://www.mathworks.com/matlabcentral/answers/348387-distinguish-uifigure-from-figure-programmatically?utm_source=chatgpt.com
y = false;
if nargin < 1, return; end
if isempty(parentfig), return; end
y = isprop(parentfig,'isUIFigure');
end
