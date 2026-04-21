function i_save2pptx(images, rmthem)


if nargin < 2, rmthem = false; end
[hasReportGen, msg] = pkg.i_isreportgenavailable('ppt');
if ~hasReportGen
    errordlg(sprintf('%s This function requires MATLAB Report Generator.', msg));
    return;
end

import mlreportgen.ppt.*;
% try
pw1 = fileparts(mfilename('fullpath'));
pth = fullfile(pw1, '..', 'assets', 'Misc', 'myTemplate.pptx');

fw = gui.gui_waitbar;
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
gui.gui_waitbar(fw);
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
