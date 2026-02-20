function i_bringtofront(fig)

if isvalid(fig)
    fig.Visible = 'on';
    fig.WindowState = 'normal';
    figure(fig);
    drawnow;
end

end
