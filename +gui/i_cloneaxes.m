    function i_cloneaxes(hAx1, hAx2)
        copyobj(allchild(hAx1), hAx2);
        
        props = {'XLim','YLim','ZLim','XScale','YScale','ZScale',...
                 'XDir','YDir','ZDir','Colormap','CLim','View',...
                 'XTicklabel','XTickLabelRotation','YTicklabel'};
        for kk = 1:numel(props)
            try
                set(hAx2, props{kk}, get(hAx1, props{kk}));
            catch
            end
        end
        
        xlabel(hAx2, get(get(hAx1,'XLabel'),'String'));
        ylabel(hAx2, get(get(hAx1,'YLabel'),'String'));
        title(hAx2,  get(get(hAx1,'Title'),'String'));
        subtitle(hAx2,  get(get(hAx1,'Subtitle'),'String'));

            gridProps = {'XGrid','YGrid','ZGrid', ...
             'XMinorGrid','YMinorGrid','ZMinorGrid', ...
             'Box'};
        for kk = 1:numel(gridProps)
            set(hAx2, gridProps{kk}, get(hAx1, gridProps{kk}));
        end

    end
