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
%   fx         - Handle to the splash screen: a gui.SplashScreen (Java),
%                a gui.NetSplashScreen (.NET), or a MATLAB figure.
%   v1         - Application version string.
%
% Three implementations, tried in order:
%   1. gui.SplashScreen, the original Java/Swing splash, when a Java
%      runtime with AWT is available.
%   2. gui.NetSplashScreen, a frameless Windows .NET splash, for Windows
%      installations without Java (MATLAB R2026b ships no JRE).
%   3. gui.sc_simplesplash, a figure-based splash, everywhere else and
%      whenever one of the above fails.

if nargin < 2, r = 0.0; end
if nargin < 1, fx = []; end

if isempty(fx)
    try
        if usejava("awt")
            [fx, v1] = in_javasplash();
        elseif ispc && NET.isNETSupported
            [fx, v1] = in_netsplash();
        else
            [fx, v1] = gui.sc_simplesplash();
        end
    catch
        [fx, v1] = gui.sc_simplesplash();
    end
else
    % Update progress on whichever splash type is active.
    v1 = '';
    if isa(fx, "gui.SplashScreen") || isa(fx, "gui.NetSplashScreen")
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
splashpng = in_pickimage();

fx = gui.SplashScreen('', char(splashpng), ...
    'ProgressBar', 'on', ...
    'ProgressPosition', 5, ...
    'ProgressRatio', 0.0);
fx.addText(30, 50, 'SCGEATOOL', 'FontSize', 18, 'Color', [1 1 1]);
fx.addText(30, 73, sprintf('Version %s', v1), ...
    'FontSize', 14, 'Color', [0.7 0.7 0.7]);
fx.addText(350, 280, 'Loading...', 'FontSize', 13, 'Color', 'white');
in_holdsplash();
end

function [fx, v1] = in_netsplash()
% Build the .NET splash, laid out like the Java one.
v1 = pkg.i_get_versionnum;
fx = gui.NetSplashScreen(in_pickimage());
fx.addText(30, 50, "SCGEATOOL", FontSize=18, Color=[1 1 1]);
fx.addText(30, 73, "Version " + v1, FontSize=14, Color=[0.7 0.7 0.7]);
fx.addText(350, 280, "Loading...", FontSize=13, Color="white");
in_holdsplash();
end

function splashpng = in_pickimage()
% Picture of the day, shared with gui.sc_simplesplash.
splashpng = i_splashimage();
if splashpng == ""
    error('gui:sc_splashscreen:noImage', 'No splash images found in assets/Images/splash_folder.');
end
end

function in_holdsplash()
% Keep the splash on screen long enough to be seen, even when the app
% initializes quickly (the caller deletes it as soon as startup finishes).
minSplashSeconds = 1.5;
drawnow;
pause(minSplashSeconds);
end
