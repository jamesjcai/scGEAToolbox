function i_bringtofront(fig)
%I_BRINGTOFRONT Show, un-minimise and raise a figure or uifigure.
%   gui.i_bringtofront(FIG) is gui.i_raisefig for a window that may be hidden
%   or minimised. Raising goes through gui.i_raisefig, because FIGURE() does
%   not raise a uifigure.

if pkg.i_isvalid(fig)
    fig.Visible = 'on';
    fig.WindowState = 'normal';
    gui.i_raisefig(fig);
    drawnow;
end

end
