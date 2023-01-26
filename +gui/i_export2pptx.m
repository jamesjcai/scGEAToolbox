function i_export2pptx(F,glist)

pw1=fileparts(mfilename('fullpath'));
pth=fullfile(pw1,'..','resources','myTemplate.pptx');
%dbfile1 = fullfile(pw1, '+run', 'external', 'stringdb', 'stringdb_human.mat');

import mlreportgen.ppt.*;
answer=questdlg('Export to PowerPoint?');
switch answer
    case 'Yes'
        fw=gui.gui_waitbar;
        if ismcc || isdeployed
            makePPTCompilable();
            warndlg('This function may not work properly in standalone applications.','');
        end
        
        N=length(F);
        images=cell(N,1);
        
        OUTppt=[tempname,'.pptx'];
        ppt = Presentation(OUTppt,pth);
        open(ppt);
        
        for k=1:N
            if isvalid(F{k})
                images{k} = [tempname,'.png'];
                saveas(F{k},images{k});
                %images{k} = [tempname,'.emf'];
                %saveas(F{k},images{k},'meta');
                
                if ~isempty(glist{1})
                    slide3 = add(ppt,'Small Title and Content');
                    replace(slide3,'Title',glist(k));
                else
                    slide3 = add(ppt,'Content Only');
                end
                replace(slide3,'Content',Picture(images{k}));
            end
        end
        close(ppt);
        rptview(ppt);
        len = length(images);
        for i = 1:len
            delete(images{i});
        end
        gui.gui_waitbar(fw);
    otherwise
end


end
