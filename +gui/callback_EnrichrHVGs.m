function callback_EnrichrHVGs(src, ~)
FigureHandle = src.Parent.Parent;
sce = guidata(FigureHandle);

gui.gui_showrefinfo('HVG Functional Analysis [PMID:31861624]');

    answer = questdlg('This function applies to a homogeneous group of cells. Continue?');
        if ~strcmp(answer, 'Yes'), return; end

            answer = questdlg('Prior to applying this function, lowly expressed genes should be filtered out and removed. Continue?');
                if ~strcmp(answer, 'Yes'), return; end

                answer = questdlg('Which HVG detecting method to use?', '', ...
                    'Brennecke et al. (2013) [PMID:24056876]', ...
                    'Splinefit Method [PMID:31697351]', ...
                    'Brennecke et al. (2013) [PMID:24056876]');
                switch answer
                    case 'Brennecke et al. (2013) [PMID:24056876]'
                        fw = gui.gui_waitbar;
                        t = sc_hvg(sce.X, sce.g, true, true);
                        gui.gui_waitbar(fw);

                        % if ~(ismcc || isdeployed)
                        %     msgfig1=export2wsdlg({'Save HVG table to variable named:'},{'T'},{t});
                        %     uiwait(msgfig1)
                        % else
                        %     gui.i_exporttable(t,true,'T');
                        % end
                        % gui.i_enrichtest(t.genes);


                    case 'Splinefit Method [PMID:31697351]'
                        fw = gui.gui_waitbar;
                        try
                            gui.sc_scatter3genes(sce.X, sce.g);
                        catch ME
                            gui.gui_waitbar(fw, true);
                            errordlg(ME.message);
                        end
                        gui.gui_waitbar(fw, true);
                    otherwise
                        return;
                end
        end
