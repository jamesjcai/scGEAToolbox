function [fx, v1] = sc_splashscreen(fx, r, ~)
% SC_SPLASHSCREEN - Display or update the application's splash screen.
%
% Usage:
%   [fx, v1] = sc_splashscreen();                % Initialize splash screen
%   sc_splashscreen(fx, r);                      % Update progress bar
%
% Inputs:
%   fx         - Handle to the splash screen (for updates).
%   r          - Progress ratio (0 to 1) for the progress bar.
%   (third arg) - Accepted but ignored; retained for signature compatibility.
%
% Outputs:
%   fx         - Handle to the splash screen (a gui.SplashScreen object when
%                Java/AWT is available, otherwise a MATLAB figure).
%   v1         - Application version string.
%
% When a Java runtime with AWT is available, this uses the original
% Java/Swing gui.SplashScreen, which renders the image and text natively and
% looks sharp. On MATLAB installations without a bundled JRE, it falls back
% to gui.sc_simplesplash, a pure-MATLAB (figure-based) splash.

if nargin < 2, r = 0.0; end
if nargin < 1, fx = []; end

if isempty(fx)
    % Initialize. Prefer the sharp Java/Swing splash when AWT is available;
    % fall back to the figure-based splash on any failure or without Java.
    if usejava("awt")
        try
            [fx, v1] = in_javasplash();
        catch
            [fx, v1] = gui.sc_simplesplash();
        end
    else
        [fx, v1] = gui.sc_simplesplash();
    end
else
    % Update progress on whichever splash type is active.
    v1 = '';
    if isa(fx, "gui.SplashScreen")
        fx.ProgressRatio = r;
    else
        gui.sc_simplesplash(fx, r);
    end
end
end

function [fx, v1] = in_javasplash()
% Build the original Java/Swing splash screen. Char literals are used for
% everything that reaches Java (image path, drawn text, color spec) because
% the underlying java.io.File / Graphics.drawString calls expect char.
v1 = pkg.i_get_versionnum;
mfolder = fileparts(mfilename('fullpath'));
splashdir = fullfile(mfolder, '..', 'assets', 'Images', 'splash_folder');

a = dir(splashdir);
a = a(~[a.isdir]); % keep files only
if isempty(a)
    error('Splash images not found in: %s', splashdir);
end

d = datetime('today');
seed = year(d) * 10000 + month(d) * 100 + day(d);
rng(seed);
idx = randi(numel(a));
splashpng = fullfile(splashdir, a(idx).name);
if ~isfile(splashpng)
    error('Splash image file not found: %s', splashpng);
end

fx = gui.SplashScreen('', splashpng, ...
    'ProgressBar', 'on', ...
    'ProgressPosition', 5, ...
    'ProgressRatio', 0.0);
fx.addText(30, 50, 'SCGEATOOL', 'FontSize', 18, 'Color', [1 1 1]);
fx.addText(30, 73, sprintf('Version %s', v1), ...
    'FontSize', 14, 'Color', [0.7 0.7 0.7]);
fx.addText(350, 280, 'Loading...', 'FontSize', 13, 'Color', 'white');

% Keep the splash on screen long enough to be seen, even when the app
% initializes quickly (the caller deletes it as soon as startup finishes).
minSplashSeconds = 1.5;
drawnow;
pause(minSplashSeconds);
end
