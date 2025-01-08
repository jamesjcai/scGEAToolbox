function [pt]=gui_3dcamera(tb, prefix, flatview)

if nargin < 3, flatview = false; end
if nargin < 2, prefix = ''; end
if nargin < 1
    hFig = gcf;
    %tb = findall(hFig, 'Type', 'uitoolbar');
    tb = uitoolbar('Parent', hFig);
    if isscalar(tb)
        tb = uitoolbar(hFig);
    else
        tb = tb(1);
    end
end
pt = uipushtool(tb, 'Separator', 'on');

try
    mfolder = fileparts(mfilename('fullpath'));
    [img, map] = imread(fullfile(mfolder, '..', 'resources', 'Images', 'camera.gif'));
    ptImage = ind2rgb(img, map);
catch
    ptImage = rand(16, 16, 3);
end
pt.CData = ptImage;
pt.Tooltip = 'Make video snapshot';
pt.ClickedCallback = @camera3dmp4;


    function camera3dmp4(~, ~)
        answer = questdlg('Make video snapshot?');
        if ~strcmp(answer, 'Yes'), return; end

        [caz,cel] = view;

        OptionZ.FrameRate = 15;
        OptionZ.Duration = 5.5;
        OptionZ.Periodic = true;
        fname = tempname;
        if ~isempty(prefix)
            [a1, b1] = fileparts(fname);
            b1 = sprintf('%s_%s', prefix, b1);
            fname = fullfile(a1, b1);
        end
        warning off
    
        if flatview
            ax = [-20, 50; -110, 65; -190, 80; -290, 60; -380, 40];
        else
            ax = [-20, 10; -110, 10; -190, 80; -290, 10; -380, 10];
        end
        CaptureFigVid(ax, fname, OptionZ);
        view(caz,cel);
        warning on
        pause(1);
        winopen(tempdir);
        pause(1);
        vfile = sprintf('%s.mp4', fname);
        if exist(vfile, 'file')
            winopen(vfile);
        else
            vfile = sprintf('%s.avi', fname);
            if exist(vfile, 'file')
                winopen(vfile);
            end
        end
        
    end
end

function CaptureFigVid(ViewZ, FileName, OptionZ)
% CaptureFigVid(ViewZ, FileName,OptionZ)
% Captures a video of the 3D plot in the current axis as it rotates based
% on ViewZ and saves it as 'FileName.mpg'. Option can be specified.
%
% ViewZ:     N-rows with 2 columns, each row are the view angles in
%            degrees, First column is azimuth (pan), Second is elevation
%            (tilt) values outside of 0-360 wrap without error,
%            *If a duration is specified, angles are used as nodes and
%            views are equally spaced between them (other interpolation
%            could be implemented, if someone feels so ambitious).
%            *If only an initial and final view is given, and no duration,
%            then the default is 100 frames.
% FileName:  Name of the file of the produced animation. Because I wrote
%            the program, I get to pick my default of mpg-4, and the file
%            extension .mpg will be appended, even if the filename includes
%            another file extension. File is saved in the working
%            directory.
% (OptionZ): Optional input to specify parameters. The ones I use are given
%            below. Feel free to add your own. Any or all fields can be
%            used
% OptionZ.FrameRate: Specify the frame rate of the final video (e.g. 30;)
% OptionZ.Duration: Specify the length of video in seconds (overrides
%    spacing of view angles) (e.g. 3.5;)
% OptionZ.Periodic: Logical to indicate if the video should be periodic.
%    Using this removed the final view so that when the video repeats the
%    initial and final view are not the same. Avoids having to find the
%    interval between view angles. (e.g. true;)
%
% % % % Example (shown in published results, video attached) % % % %
% figure(171);clf;
% surf(peaks,'EdgeColor','none','FaceColor','interp','FaceLighting','phong')
% daspect([1,1,.3]);axis tight;
% OptionZ.FrameRate=15;OptionZ.Duration=5.5;OptionZ.Periodic=true;
% CaptureFigVid([-20,10;-110,10;-190,80;-290,10;-380,10],'WellMadeVid',OptionZ)
%
% Known issues: MPEG-4 video option only available on Windows machines. See
% fix where the VideoWriter is called.
%
% Getframe is used to capture image and current figure must be on monitor 1
% if multiple displays are used. Does not work if no display is used.
%
% Active windows that overlay the figure are captured in the movie.  Set up
% the current figure prior to calling the function. If you don't specify
% properties, such as tick marks and aspect ratios, they will likely change
% with the rotation for an undesirable effect.

% Cheers, Dr. Alan Jennings, Research assistant professor,
% Department of Aeronautics and Astronautics, Air Force Institute of Technology
if nargin < 3
    OptionZ = struct([]);
end

% check orientation of ViewZ, should be two columns and >=2 rows
if size(ViewZ, 2) > size(ViewZ, 1)
    ViewZ = ViewZ.';
end
if size(ViewZ, 2) > 2
    warning('AJennings:VidWrite', ...
        'Views should have n rows and only 2 columns. Deleting extraneous input.');
    ViewZ = ViewZ(:, 1:2); %remove any extra columns
end

if ispc
    daObj = VideoWriter(FileName, 'MPEG-4'); %my preferred format
else
    daObj=VideoWriter(FileName); %for default video format.
end
% MPEG-4 CANNOT BE USED ON UNIX MACHINES
% set values:
% Frame rate
if isfield(OptionZ, 'FrameRate')
    daObj.FrameRate = OptionZ.FrameRate;
end
if isfield(OptionZ, 'Duration') %space out view angles
    temp_n = round(OptionZ.Duration*daObj.FrameRate); % number frames
    temp_p = (temp_n - 1) / (size(ViewZ, 1) - 1); % length of each interval
    ViewZ_new = zeros(temp_n, 2);
    for inis = 1:(size(ViewZ, 1) - 1)
        ViewZ_new(round(temp_p*(inis - 1)+1):round(temp_p*inis+1), :) = ...
            [linspace(ViewZ(inis, 1), ViewZ(inis+1, 1), ...
            round(temp_p*inis)-round(temp_p*(inis - 1))+1).', ...
            linspace(ViewZ(inis, 2), ViewZ(inis+1, 2), ...
            round(temp_p*inis)-round(temp_p*(inis - 1))+1).'];
    end
    ViewZ = ViewZ_new;
end
if length(ViewZ) == 2 % only initial and final given
    ViewZ = [linspace(ViewZ(1, 1), ViewZ(end, 1)).', ...
        linspace(ViewZ(1, 2), ViewZ(end, 2)).'];
end
if isfield(OptionZ, 'Periodic') && OptionZ.Periodic == true
    ViewZ = ViewZ(1:(end -1), :); %remove last sample
end
open(daObj);
for kathy = 1:size(ViewZ, 1)
    view(ViewZ(kathy, :));
    drawnow;
    writeVideo(daObj, getframe(gcf)); %use figure, since axis changes size based on view
end
close(daObj);
end