function [fx] = sc_splashscreen(fx, r, showloading)

if nargin<3, showloading=true; end
if nargin<2, r = 0.0; end
if nargin<1
    mfolder = fileparts(mfilename('fullpath'));
    v1='v24.3.3';
    %[~, v1] = pkg.i_majvercheck;
    splashpng = '700813831-hero-1536x1536.png';
    %splashpng = 'Untitled.png';
    fx = gui.SplashScreen( 'Splashscreen', ...
        fullfile(mfolder,'..','resources', splashpng), ...
                    'ProgressBar', 'on', ...
                    'ProgressPosition', 5, ...
                    'ProgressRatio', 0.0 );
    fx.addText( 30, 50, 'SCGEATOOL', 'FontSize', 18, 'Color', [1 1 1] )
    fx.addText( 30, 80, v1, 'FontSize', 14, 'Color', [0.7 0.7 0.7] )
    if showloading
        fx.addText( 307, 277, 'Loading...', 'FontSize', 14, 'Color', 'white' )
    end
end
if nargin==2
    fx.ProgressRatio=r;
end