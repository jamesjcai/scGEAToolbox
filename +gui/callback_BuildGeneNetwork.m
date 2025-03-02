function callback_BuildGeneNetwork(src, ~)
    [FigureHandle, sce] = gui.gui_getfigsce(src);
    
    [glist] = gui.i_selectngenes(sce, [], FigureHandle);
    if isempty(glist), return; end
    
    [y, i] = ismember(upper(glist), upper(sce.g));
    if ~all(y), error('xxx'); end
    fprintf("%s\n", glist)
    
    [Xt] = gui.i_transformx(sce.X, true, 5);
    if isempty(Xt), return; end
    
    x = Xt(i, :);

    switch questdlg("Select algorithm:",'',"PC Regression","Chaterjee Correlation","PC Regression")
        case "PC Regression"
            fw = gui.gui_waitbar;
            A = sc_pcnet(x);
            gui.gui_waitbar(fw);
        case "Chaterjee Correlation"
            fw = gui.gui_waitbar_adv;
            n = size(x,1);
            A = zeros(n);
            for k=1:n
                gui.gui_waitbar_adv(fw,k/n);
                for l=1:n
                    if k~=l
                        A(k,l) = pkg.e_xicor(x(k,:),x(l,:));
                    end
                end
            end
            gui.gui_waitbar_adv(fw);
        otherwise
            return;
    end
%    [~, systemView] = memory;
% disp(systemView.PhysicalMemory.Available)
%    bytesPerElement = 8;    % For double precision
%    maxElements = systemView.PhysicalMemory.Available / bytesPerElement;
%    maxSize = floor(sqrt(maxElements));  % If square matrix is needed


    sc_grnview(A, glist, [], FigureHandle);
end
