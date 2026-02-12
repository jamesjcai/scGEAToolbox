function callback_PickColorMap(src, n, showzero, showdlg)

if isa(src, 'matlab.ui.Figure')
    FigureHandle = src;
else
    [FigureHandle, sce] = gui.gui_getfigsce(src);
    n = numel(unique(sce.c));
end

    

    if nargin < 4 || isempty(showdlg), showdlg = false; end 
    if nargin < 3 || isempty(showzero), showzero = false; end
    
    
    % disp(sprintf('Using %d colors',n));
    % n = max([n, 3]);
    folder = fileparts(mfilename('fullpath'));
    %a = strfind(folder, filesep);
    %folder = extractBefore(folder, a(end)+1);
    wrkpth = fullfile(folder, '..', 'external', 'ml_cbrewer');

    if ~(ismcc || isdeployed)
        addpath(wrkpth);
    end
    CT = cbrewer('seq', 'Blues', max([n, 3]));
    
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
        linspecer(n, 'red'), linspecer(n, 'gray'), ...
        linspecer(n, 'green'), ...
        purpleblue(n), ...
        purplemap(n),...
        seuratFeaturePlot(n)};
    cn = {'autumnzero', '7lines', 'parula', 'summer', 'jet', 'copper', ...
        'winter', 'hsv', ...
        'UKRAINE', ...
        'mycmap', 'seqBlues', ...
        'divSpectral', 'divRdBu', 'seqPuBuGn', 'qualSet1', 'qualDark2', ...
        '12qualLinspecer', 'seqLinespecer', 'redLinespecer', 'grayLinespecer', ...
        'greenLinespecer',...
        'purpleblue','purplemap','seuratFeaturePlot'};
    assert(length(co) == length(cn));
    


    if showdlg
        [indx, tf] = gui.myListdlg(FigureHandle, cn, 'Pick a colormap', ...
                    [], false);
        if ~tf
            retrun;
        end
    else
        indx = randi(length(co));
    end
    
    colormap(FigureHandle, abs(co{indx}));
    % fprintf('Set colormap to %s.\n', cn{indx});
    
    if showzero
        cm = colormap;
        cm(1, :) = [.8, .8, .8];
        colormap(FigureHandle, cm);
    end
end



function cmap = purpleblue(n)

if nargin < 1
    n = 256;
end

anchors = [
    0.92 0.88 1.00
    0.75 0.70 0.98
    0.55 0.55 0.98
    0.35 0.40 0.95
    0.10 0.20 0.90
];

x  = linspace(0,1,size(anchors,1));
xi = linspace(0,1,n);

cmap = interp1(x, anchors, xi, 'pchip');
end

function cmap = purplemap(n)
%PURPLEMAP  Lavender -> deep purple colormap

if nargin < 1
    n = 256;
end

anchors = [
    0.93 0.88 0.98   % very light lavender
    0.80 0.70 0.95
    0.62 0.55 0.96
    0.42 0.38 0.94
    0.22 0.18 0.80   % deep purple
];

x  = linspace(0,1,size(anchors,1));
xi = linspace(0,1,n);

cmap = interp1(x, anchors, xi, 'pchip');
end

function cmap = seuratFeaturePlot(n)

if nargin < 1
    n = 256;
end

anchors = [
    0.83 0.83 0.83   % lightgrey  (Seurat low)
    0.00 0.00 1.00   % blue       (Seurat high)
];

x  = linspace(0,1,size(anchors,1));
xi = linspace(0,1,n);

cmap = interp1(x, anchors, xi, 'pchip');

end


