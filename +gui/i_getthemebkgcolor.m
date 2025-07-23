function bgColor = i_getthemebkgcolor(fig)
    arguments
        fig = gcf % Default to current figure
    end
    bgColor = [1 1 1];
    
    try
        % Ensure the figure has theme support
        gt = theme(fig); % Returns a GraphicsTheme object :contentReference[oaicite:1]{index=1}
    
        % Determine if the theme is light or dark
        style = gt.BaseColorStyle;  % "light" or "dark" :contentReference[oaicite:2]{index=2}
    
        % Ask MATLAB for the actual background color
        % bgColor = fig.Color;  % RGB triplet
        if strcmp(style, 'dark')
            bgColor = [0 0 0]; % Set to white for light theme
        end
    catch
    end
end
