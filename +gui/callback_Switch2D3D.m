function []=callback_Switch2D3D(src, ~)      % xxx
    [para] = gui.i_getoldsettings(src);
    if isempty(h.ZData)               % current 2D
        if ~(size(sce.s, 2) >= 3)
            in_EmbeddingAgain(src, [], 3);
        else
            h = gui.i_gscatter3(sce.s, c, methodid, hAx);
            if ~isempty(ax) && ~isempty(bx) && ~any([ax, bx] == 0)
                view(ax, bx);
            else
                view(3);
            end
        end
    else        % current 3D do following
        answer = questdlg('Select a method of making 2D embedding.','', ...
            'Re-embedding to 2D','Compress 3D embedding','Cancel', ...
            'Re-embedding to 2D');
        switch answer
            case 'Cancel'
                return;
            case 'Re-embedding to 2D'
                in_EmbeddingAgain(src, [], 2);
            case 'Compress 3D embedding'
                [ax, bx] = view();
                answer2 = questdlg('Which view to be used to project cells?', '', ...
                    'X-Y Plane', 'Screen/Camera', 'PCA-rotated', 'X-Y Plane');
                switch answer2
                    case 'X-Y Plane'
                        sx = sce.s;
                    case 'Screen/Camera'
                        sx = pkg.i_3d2d(sce.s, ax, bx);
                    case {'PCA-rotated'}
                        [~, sx] = pca(sce.s);
                    otherwise
                        return;
                end
                h = gui.i_gscatter3(sx(:, 1:2), c, methodid, hAx);
                sce.s = sx;
        end
    end
    title(sce.title);
    subtitle('[genes x cells]');

    h.Marker = para.oldMarker;
    h.SizeData = para.oldSizeData;
    colormap(para.oldColorMap);
end
