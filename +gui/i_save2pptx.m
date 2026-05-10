function i_save2pptx(images, rmthem, fw, parentfig)


if nargin < 2, rmthem = false; end
if nargin < 3, fw = []; end
if nargin < 4, parentfig = []; end
[hasReportGen, msg] = pkg.i_isreportgenavailable('ppt');
if ~hasReportGen
    errordlg(sprintf('%s This function requires MATLAB Report Generator.', msg));
    return;
end

import mlreportgen.ppt.*;
% try
pw1 = fileparts(mfilename('fullpath'));
pth = fullfile(pw1, '..', 'assets', 'Misc', 'myTemplate.pptx');

ownsWaitbar = nargin < 3 || isempty(fw);
if ownsWaitbar
    fw = gui.gui_waitbar;
else
    gui.myWaitbar(parentfig, fw, false, '', 'Exporting PowerPoint...', 0.995);
end
OUTppt = [tempname, '.pptx'];
ppt = Presentation(OUTppt, pth);
open(ppt);
for k = 1:length(images)
        slide3 = add(ppt, 'Content Only');
        % slide3 = add(ppt,'Small Title and Content');
        % replace(slide3,'Title',glist(k));
        replace(slide3, 'Content', Picture(images{k}));
    end
    % pictureSlide = add(ppt,'Title and Picture',2);
close(ppt);
if ownsWaitbar
    gui.gui_waitbar(fw);
end
pkg.i_openoutputfile(OUTppt);
% catch ME
%     gui.gui_waitbar(fw, true);
%     errordlg(ME.message);
% end

if rmthem
    len = length(images);
    for i = 1:len
        delete(images{i});
    end
end
