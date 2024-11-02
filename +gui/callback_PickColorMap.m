function callback_PickColorMap(src, ~, showzero)

FigureHandle = src.Parent.Parent;
sce = guidata(FigureHandle);
if nargin < 3, showzero = false; end
n = numel(unique(sce.c));
% disp(sprintf('Using %d colors',n));
n = max([n, 3]);
folder = fileparts(mfilename('fullpath'));
a = strfind(folder, filesep);
folder = extractBefore(folder, a(end)+1);
wrkpth = fullfile(folder, '+run', 'external', 'ml_cbrewer');
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


%     co={cx,a,b,'default',summer(kc),...
%         jet(kc),copper(kc),winter(kc)};

indx = randi(length(co));
colormap(abs(co{indx}));
fprintf('Set colormap to %s.\n', cn{indx});

% answer=questdlg('Random colormap?','');
% switch answer
%     case 'Yes'
%         indx=randi(length(co));
%         colormap(co{indx});
%         % mb=helpdlg(sprintf('Using colormap: %s',cn{indx}),'');
%
%     case 'No'
%         [indx,tf] = listdlg('PromptString',{'Select a colormap'},...
%             'SelectionMode','single','ListString',cn);
%         if tf==1
%             colormap(co{indx});
%         end
%     case 'Cancel'
%         return;
%     otherwise
%         return;
% end

if showzero
    cm = colormap;
    cm(1, :) = [.8, .8, .8];
    colormap(cm);
end


end
