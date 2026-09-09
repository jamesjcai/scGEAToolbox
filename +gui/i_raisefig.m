function i_raisefig(parentfig)
% i_raisefig - Bring a figure or uifigure back to the front.
%
% i_raisefig(PARENTFIG) raises PARENTFIG above other windows. uifigures are
% raised with focus(), which is the supported way for App Designer windows;
% figure() does not even make a uifigure the groot CurrentFigure. Classic
% figures are raised with figure(). An empty or deleted handle is ignored, so
% this is safe to call from onCleanup blocks and dialog teardown code.

if nargin < 1 || isempty(parentfig) || ~pkg.i_isvalid(parentfig)
    return;
end

if gui.i_isuifig(parentfig)
    try
        focus(parentfig);
    catch
        % focus() is unavailable on older releases; nothing else to try
    end
else
    figure(parentfig);
end
end
