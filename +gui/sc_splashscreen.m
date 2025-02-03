function [fx, v1] = sc_splashscreen(fx, r, showloading)
% SC_SPLASHSCREEN - Display or update the application's splash screen.
% 
% Usage:
%   [fx, v1] = sc_splashscreen();                % Initialize splash screen
%   sc_splashscreen(fx, r);                      % Update progress bar
%
% Inputs:
%   fx         - Handle to the splash screen object (for updates).
%   r          - Progress ratio (0 to 1) for the progress bar.
%   showloading - (Optional) Boolean to show "Loading..." text (default: true).
%
% Outputs:
%   fx         - Handle to the splash screen object.
%   v1         - Application version string.
% ```

import pkg.*
import gui.*
if nargin<3, showloading = true; end
if nargin<2, r = 0.0; end
if nargin<1
    mfolder = fileparts(mfilename('fullpath'));
    v1 = pkg.i_get_versionnum;
    pngfilename = 'splash.png';
    splashpng = fullfile(mfolder, '..','resources', 'Images', pngfilename);
    if ~isfile(splashpng)
        error('Splash image file not found: %s', splashpng);
    end    
    fx = gui.SplashScreen('', splashpng, ...
                    'ProgressBar', 'on', ...
                    'ProgressPosition', 5, ...
                    'ProgressRatio', 0.0 );
    textX = 30; textY = 50;
    fx.addText(textX, textY, 'SCGEATOOL', 'FontSize', 18, 'Color', [1 1 1]);
    fx.addText(30, 73, sprintf('Version %s', v1), ...
        'FontSize', 14, 'Color', [0.7 0.7 0.7] )
    if showloading
        fx.addText(350, 280, 'Loading...', ...
            'FontSize', 13, 'Color', 'white')
    end
end
if nargin == 2
    fx.ProgressRatio = r;
end