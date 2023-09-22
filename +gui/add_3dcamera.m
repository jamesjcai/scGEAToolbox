function add_3dcamera(tb, prefix, flatview)

if nargin < 3, flatview = false; end
if nargin < 2, prefix = ''; end
if nargin < 1
    hFig = gcf;
    tb = findall(hFig, 'Type', 'uitoolbar');
    if length(tb) == 1
        tb = uitoolbar(hFig);
    else
        tb = tb(1);
    end
end
pt = uipushtool(tb, 'Separator', 'on');
folder = fileparts(mfilename('fullpath'));
a = strfind(folder, filesep);
folder = extractBefore(folder, a(end)+1);
wrkpth = fullfile(folder, 'resources', 'camera.gif');
[img, map] = imread(wrkpth);
ptImage = ind2rgb(img, map);
pt.CData = ptImage;
pt.Tooltip = 'Make video snapshot';
pt.ClickedCallback = @camera3dmp4;


    function camera3dmp4(~, ~)
        answer = questdlg('Make video snapshot?');
        if ~strcmp(answer, 'Yes'), return; end
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
            gui.CaptureFigVid([-20, 50; -110, 65; -190, 80; -290, 60; -380, 40], fname, OptionZ);
        else
            gui.CaptureFigVid([-20, 10; -110, 10; -190, 80; -290, 10; -380, 10], fname, OptionZ);
        end

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