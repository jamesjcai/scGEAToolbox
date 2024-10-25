function [fx, v1] = sc_splashscreen(fx, r, showloading)
import pkg.*
import gui.*
if nargin<3, showloading=true; end
if nargin<2, r = 0.0; end
if nargin<1
    mfolder = fileparts(mfilename('fullpath'));
    v1 = pkg.i_get_versionnum;
    %pngfilename = 'dna-adn-black-background-1080P-wallpaper.jpg';
    %pngfilename = 'OIP (1).jpg';
    %pngfilename = 'wave-white-particles-abstract-technology-flow-background-future-vector-illustration_435055-172.jpg';
    %pngfilename = '0321efeecd4da13290687593d55674cf.png';
    pngfilename = '700813831-hero-1536x1536.png';
    splashpng = fullfile(mfolder,'..','resources', 'Images', pngfilename);
    %splashpng = '09-small.jpg';
    %splashpng = 'Untitled.png';
    fx = gui.SplashScreen('', splashpng, ...
                    'ProgressBar', 'on', ...
                    'ProgressPosition', 5, ...
                    'ProgressRatio', 0.0 );
    fx.addText( 30, 50, 'SCGEATOOL', 'FontSize', 18, 'Color', [1 1 1] )
    fx.addText( 30, 73, sprintf('Version %s', v1), ...
        'FontSize', 14, 'Color', [0.7 0.7 0.7] )
    if showloading
        fx.addText( 307, 277, 'Loading...', 'FontSize', 14, 'Color', 'white' )
    end
end
if nargin==2
    fx.ProgressRatio=r;
end