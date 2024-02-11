function pt=ui2getpushbuttoncddate
mfolder = fileparts(mfilename('fullpath'));

FigureHandle = uifigure;
MainToolbarHandle = uitoolbar(FigureHandle);

        pt = uipushtool(MainToolbarHandle, 'Separator', 'on');
        imgFil='plotpicker-pzmap.gif';
        pt.Icon = fullfile('..','resources',imgFil);
        %pt.CData = in_getPtImage(imgFil);
        pt.Tooltip = 'tooltipTxt';
        %pt.CData

        pt.Icon

    function [ptImage] = in_getPtImage(imgFil)
        try
            [img, map] = imread(fullfile(mfolder,'..', 'resources', imgFil));
            ptImage = ind2rgb(img, map);
        catch
            try
                [img, map] = imread(fullfile(matlabroot,'toolbox', ...
                    'matlab','icons', imgFil));
                 ptImage = ind2rgb(img, map);
            catch
                ptImage = rand(16, 16, 3);
            end
        end
    end
end
        

