function varargout = animsplash(varargin)
%ANIMSPLASH  Creates a gif splash screen.
%   HSPLASH =  ANIMSPLASH(FILENAME,HEIGHT,WIDTH) creates a gif splash screen and returns its handle
%   in HSPLASH.Both HEIGHT and WIDTH are pixelssize of GIF. The gif splash screen will be always shown 
%   until ANIMSPLASH(HSPLASH,'off') is called to turn it off.
%
%   ANIMSPLASH(HSPLASH,'off') closes the gif splash screen with handle HSPLASH. Since
%   HSPLASH is actually a javax.swing.JWindow object, you can also close it
%   by HSPLASH.dispose();
%
%
%   Example
%       % To show a gif splash screen while you GUI program is initialized
%       h = animsplash('anim.gif',400,200);
%                       ^could also be 'Data/anim.gif' if GIF is in a Folder
%                       400 = Height of GIF in pixels
%                       200 = Width of GIF in pixels
%       ............;% Place you GUI initializing code here
%       animsplash(h,'off'); % Close the gif splash screen
%       ............; % Other commands


%   Stefan Braun, April 2016
%   Revision: 1.0  Date: 26/04/2016

%% Check output
if nargout >=2  
    error('MATLAB:animsplash','%s','Too many output!');
end
%% Check input
[filename, Height, Width, handle, msg] = parse_inputs(varargin{:});
if (~isempty(msg))
    error('MATLAB:splash:inputParsing', '%s', msg);
end

if (~isempty(handle))
    handle.dispose;
    return;
end


%% Create splash screen
win = javax.swing.JWindow;
html=strcat('<html><img src="file:./',filename,'"/></html>');
label = javax.swing.JLabel(html);
win.getContentPane.add(label);
win.setAlwaysOnTop(true);
win.pack;

%% set the splash image to the center of the screen
screenSize = win.getToolkit.getScreenSize;
screenHeight = screenSize.height;
screenWidth = screenSize.width;
% set Location
win.setLocation((screenWidth-Width)/2,(screenHeight-Height)/2);

win.show % show the splash screen

%% Output the handle
if (nargout==1)
    varargout{1} = win;
end



%% Function parse_inputs
function [filename, Height, Width, handle, msg] = parse_inputs(varargin)

filename = '';
handle = [];
Height = 400;
Width = 400;
msg = '';

% Parse arguments based on their number.
switch(nargin)
    case 0  % Not allowed.
        msg = 'Too few input arguments.';
        return;
    case 1  % Filename.
        msg = 'Too few input arguments.';
    case 2  % Filename + handle 'off'
        in1 = varargin{1};
        in2 = varargin{2};
        if isjava(in1) && isequal(in2,'off')
            handle = in1;
        else
            msg='Input type mismatch. Help animsplash for more information';
        end
    case 3 
        in1 = varargin{1};
        in2 = varargin{2};
        in3 = varargin{3};
        if ischar(in1) && isnumeric(in2) && isnumeric(in3)
            filename = in1;
            Height = in2;
            Width = in3;
        else
            msg='Input type mismatch. Help animsplash for more information';
        end
    otherwise
        msg = 'Too many input arguments.';
end