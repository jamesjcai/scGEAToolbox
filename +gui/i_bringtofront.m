function i_bringtofront(fig)

if pkg.i_isvalid(fig)
    fig.Visible = 'on';
    fig.WindowState = 'normal';
    figure(fig);
    drawnow;
end

end
