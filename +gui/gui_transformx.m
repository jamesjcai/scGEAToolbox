function [X]=gui_transformx(X)
    if nargin<1
        X=nbinrnd(20,0.98,1000,200);
        disp('Using simulated X.');
    end
    
    answer = questdlg('Normalize, transform or impute X?');
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
    
    listitems={'(1) Library Size Normalization',...
        '(2) Log10(x+1) Transformation',...
        '(1)+(2)',...
        'Pearson Residual Transformation',...
        'MAGIC Imputation'};
    [indx,tf] = listdlg('PromptString',{'Select Method',...
        '',''},'SelectionMode','single','ListString',listitems,'ListSize',[200 300]);
    if tf==1
        fw=gui.gui_waitbar;
        try
        switch indx
            case 1
                X=sc_norm(X);
            case 2
                X=log10(X+1);
            case 3
                X=sc_norm(X);
                X=log10(X+1);
            case 4
                X=sc_transform(X);
            case 5
                X=run.MAGIC(X,true);
        end
        catch ME
            gui.gui_waitbar(fw);
            errordlg(ME.message)
            rethrow(ME)
        end
        gui.gui_waitbar(fw);
    end
end

