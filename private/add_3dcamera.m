function add_3dcamera(tb)
if nargin<1
hFig=gcf;
tb = findall(hFig,'Type','uitoolbar');
if length(tb)==1
    tb = uitoolbar(hFig);
else
    tb=tb(1);
end

end
pt = uipushtool(tb,'Separator','on');
[img,map] = imread(fullfile(fileparts(which(mfilename)),...
            'camera.gif'));
ptImage = ind2rgb(img,map);
pt.CData = ptImage;
pt.Tooltip = 'Select a gene to show expression';
pt.ClickedCallback = @camera3dmp4;


function camera3dmp4(~,~)
    OptionZ.FrameRate=15;OptionZ.Duration=5.5;OptionZ.Periodic=true;
    fname=tempname;
    CaptureFigVid([-20,10;-110,10;-190,80;-290,10;-380,10],fname,OptionZ);
    pause(1);
    winopen(tempdir);
    pause(3);
    winopen(sprintf('%s.mp4',fname));
end

end