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
pw1 = fileparts(mfilename('fullpath'));
pth = fullfile(pw1, '..', 'assets', 'Misc', 'myTemplate.pptx');

ownsWaitbar = nargin < 3 || isempty(fw);
if ownsWaitbar
    fw = gui.myWaitbar([]);
else
    gui.myWaitbar(parentfig, fw, false, '', 'Exporting PowerPoint...', 0.995);
end
% Same path as GUI.I_EXPORT2PPTX: a named per-process folder and a
% timestamped file, so the deck can be named in the fallback message
% below when the shell declines to open it.
[~, OUTppt] = pkg.i_tempdirfile('scgeatool_pptx', 'pptx');
ppt = Presentation(OUTppt, pth);
open(ppt);
for k = 1:length(images)
        slide3 = add(ppt, 'Content Only');
        replace(slide3, 'Content', Picture(images{k}));
    end
close(ppt);
if ownsWaitbar
    gui.myWaitbar([], fw);
end
if ~pkg.i_openoutputfile(OUTppt)
    gui.myHelpdlg(parentfig, sprintf( ...
        ['The presentation could not be opened here. ', ...
        'It has been saved as\n\n%s'], OUTppt));
end

if rmthem
    len = length(images);
    for i = 1:len
        delete(images{i});
    end
end
