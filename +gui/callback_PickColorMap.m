function callback_PickColorMap(src, ~, showzero)

    [FigureHandle, sce] = gui.gui_getfigsce(src);

    if nargin < 3, showzero = false; end
    n = numel(unique(sce.c));
    % disp(sprintf('Using %d colors',n));
    n = max([n, 3]);
    folder = fileparts(mfilename('fullpath'));
    %a = strfind(folder, filesep);
    %folder = extractBefore(folder, a(end)+1);
    wrkpth = fullfile(folder, '..', 'external', 'ml_cbrewer');

    wrkpth

    if ~(ismcc || isdeployed)
        addpath(wrkpth);
    end
    CT = cbrewer('seq', 'Blues', n);
    
    cx = autumn(n);
    cx(1, :) = [.8, .8, .8];
    %a=lines(kc);
    %rng("shuffle");
    %b=a(randperm(size(a,1)),:);
    ukraineflag = [0, 87, 183; 255, 215, 0] ./ 255;
    mycmap = pkg.i_mycolormap(n);
    co = {cx, lines(n), parula(n), summer(n), ...
        jet(n), copper(n), winter(n), hsv(n), ...
        ukraineflag, ...
        mycmap, CT, ...
        cbrewer('div', 'Spectral', n), ...
        cbrewer('div', 'RdBu', n), ...
        cbrewer('seq', 'PuBuGn', n), ...
        cbrewer('qual', 'Set1', n), ...
        cbrewer('qual', 'Dark2', n), ...
        linspecer(min([n, 12]), 'qualitative'), ...
        linspecer(n, 'sequential'), ...
        linspecer(n, 'red'), linspecer(n, 'gray'), linspecer(n, 'green')};
    cn = {'autumnzero', '7lines', 'parula', 'summer', 'jet', 'copper', ...
        'winter', 'hsv', ...
        'UKRAINE', ...
        'mycmap', 'seqBlues', ...
        'divSpectral', 'divRdBu', 'seqPuBuGn', 'qualSet1', 'qualDark2', ...
        '12qualLinspecer', 'seqLinespecer', 'redLinespecer', 'grayLinespecer', ...
        'greenLinespecer'};
    assert(length(co) == length(cn));
    
    indx = randi(length(co));
    colormap(FigureHandle, abs(co{indx}));
    fprintf('Set colormap to %s.\n', cn{indx});
    
    if showzero
        cm = colormap;
        cm(1, :) = [.8, .8, .8];
        colormap(FigureHandle, cm);
    end
end