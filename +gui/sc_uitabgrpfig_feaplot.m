function sc_uitabgrpfig_feaplot(feays, fealabels, sce_s, ...
                parentfig, methodid, cazcel)

if nargin < 6, cazcel = []; end
if nargin < 5, methodid = 1; end
if nargin < 4, parentfig = []; end

if ~isstring(fealabels), fealabels = string(fealabels); end

if (ismcc || isdeployed) && pkg.i_isreportgenavailable('ppt'), makePPTCompilable(); end


hx = gui.myFigure(parentfig);
hFig = hx.FigHandle;

n = length(fealabels);

tabgp = uitabgroup();
tab = cell(n,1);
ax0 = cell(n,1);
ax = cell(n,2);

idx = 1;
focalg = fealabels(idx);

for k=1:n
    c = feays{k};
    if issparse(c), c = full(c); end
    if ~isnumeric(c)
        [c] = findgroups(string(c));
    end
    tab{k} = uitab(tabgp, 'Title', sprintf('%s', fealabels(k)));

    ax0{k} = axes('parent', tab{k});
    ax{k,1} = ax0{k};

    switch methodid
        case 1
            if size(sce_s,2) > 2
                scatter3(sce_s(:,1), sce_s(:,2), sce_s(:,3), 5, c, 'filled');
            else
                scatter(sce_s(:,1), sce_s(:,2), 5, c, 'filled');
            end
            if ~isempty(cazcel)
                view(ax{k,1}, [cazcel(1), cazcel(2)]);
            end
        case 2
            gui.i_stemscatter(sce_s, feays{k});
            title(ax{k,1}, strrep(fealabels(k),'_','\_'));
    end


end

tabgp.SelectionChangedFcn=@displaySelection;

hx.addCustomButton('off',  @i_genecards, 'www.jpg', 'GeneCards...');
hx.addCustomButton('off', @i_proteinstructure, 'hexagon_16dp_000000_FILL0_wght400_GRAD0_opsz20.jpg', 'Protein Structure...');

hx.addCustomButton('off', @in_savedata, "floppy-disk-arrow-in.jpg", 'Save Gene List...');
hx.show(parentfig);


function in_savedata(~,~)
        gui.i_exporttable(table(fealabels), true, ...
            'Tmarkerlist', 'MarkerListTable', [], [], hFig);
    end

function displaySelection(~,event)
        t = event.NewValue;
        txt = t.Title;
        [~,idx]=ismember(txt,fealabels);
        focalg = fealabels(idx);
    end

function i_genecards(~, ~)
        web(sprintf('https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s', focalg),'-new');
    end

function i_proteinstructure(~, ~)
        gui.i_viewprotein(focalg, ParentFig=hFig);
    end
end
