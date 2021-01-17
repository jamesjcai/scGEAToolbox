function [X]=gui_transformx(X)
    if nargin<1
        X=nbinrnd(20,0.98,1000,200);
        disp('Using simulated X.');
    end
    
    answer = questdlg('Transform X?');
    if strcmp(answer,'Yes')
        % 
    elseif strcmp(answer,'No')
        return;
    elseif strcmp(answer,'Cancel')
        X=[];
        return;
    else
        error('Wrong option');        
    end
    
    listitems={'Library Size','Pearson Residuals',...
        'MAGIC'};
    [indx,tf] = listdlg('PromptString',{'Select Transformation Method',...
    '',''},'SelectionMode','single','ListString',listitems);
    if tf==1        
        switch indx
            case 1
                X=sc_norm(X);
            case 2
                X=sc_transform(X);
            case 3
                X=run.MAGIC(X,true);
        end
    else
        return;
    end
end

