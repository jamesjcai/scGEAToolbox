function i_pickcolormap(~,~,c,setzerocolor)
    if nargin<4, setzerocolor=false; end
    if setzerocolor
        list = {'parula','turbo','hsv','hot','cool','spring',...
                'summer','autumn (default)',...
                'winter','jet'};
    else
        list = {'parula','turbo','hsv','hot','cool','spring',...
                'summer','autumn',...
                'winter','jet','gray','bone','pink','copper'};
    end
    [indx,tf] = listdlg('ListString',list,'SelectionMode','single',...
                        'PromptString','Select a colormap:');
    if tf==1
        a=list{indx};
        if strcmp(a,'autumn (default)')
            a='autumn';
        end
        if setzerocolor
            gui.i_setautumncolor(c,a);
        else
            colormap(a);
        end
        % setpref('scgeatoolbox','prefcolormapname',a);        
    end
end