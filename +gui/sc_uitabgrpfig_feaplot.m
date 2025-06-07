function sc_uitabgrpfig_feaplot(feays, fealabels, sce_s, ...
                parentfig, methodid, cazcel)

if nargin < 6, cazcel = []; end
if nargin < 5, methodid = 1; end
if nargin < 4, parentfig = []; end

if ~isstring(fealabels), fealabels = string(fealabels); end

if ismcc || isdeployed, makePPTCompilable(); end
import mlreportgen.ppt.*;

pw1 = fileparts(mfilename('fullpath'));
pth = fullfile(pw1, '..', 'assets', 'Misc', 'myTemplate.pptx');


hx=gui.myFigure;
hFig = hx.FigHandle;
% hFig.Position(3) = hFig.Position(3) * 1.8;

n = length(fealabels);
%a = getpref('scgeatoolbox', 'prefcolormapname', 'autumn');

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
        [c] = grp2idx(c);
    end
    tab{k} = uitab(tabgp, 'Title', sprintf('%s', fealabels(k)));
    
    ax0{k} = axes('parent', tab{k});
    %ax{k,1} = subplot(1,2,1);
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
            zlabel(strrep(fealabels(k),'_','\_'));
    end


    % ax{k,2} = subplot(1,2,2);
    % scatter(sce_s(:,1), sce_s(:,2), 5, c, 'filled');
    % stem3(sce_s(:,1), sce_s(:,2), c, 'marker', 'none', 'color', 'm');
    % hold on;
    % scatter3(sce_s(:,1), sce_s(:,2), zeros(size(sce_s(:,2))), 5, c, 'filled');
    % title(ax{k,1}, strrep(fealabels(k),'_','\_'));
    % subtitle(ax{k,1}, gui.i_getsubtitle(c));
    % title(ax{k,2}, strrep(fealabels(k),'_','\_'));
    % subtitle(ax{k,2}, gui.i_getsubtitle(c));
    % gui.i_setautumncolor(c, a, true, any(c==0));
end
  
tabgp.SelectionChangedFcn=@displaySelection;

hx.addCustomButton('off',  @i_genecards, 'www.jpg', 'GeneCards...');
%hx.addCustomButton('on', {@i_PickColorMap, c}, 'plotpicker-compass.gif', 'Pick new color map...');
%hx.addCustomButton('off', @i_RescaleExpr, 'IMG00074.GIF', 'Rescale expression level [log2(x+1)]');
%hx.addCustomButton('off', @i_ResetExpr, 'plotpicker-geobubble2.gif', 'Reset expression level');
%hx.addCustomButton('off', {@gui.i_savemainfig, 3}, "powerpoint.gif", 'Save Figure to PowerPoint File...');

hx.addCustomButton('off', @in_savedata, "floppy-disk-arrow-in.jpg", 'Save Gene List...');
hx.show(parentfig);


    function in_savedata(~,~)
        gui.i_exporttable(table(fealabels), true, ...
            'Tmarkerlist','MarkerListTable');    
    end
    
    function displaySelection(~,event)
        t = event.NewValue;
        txt = t.Title;
        % disp("Viewing gene " + txt);
        [~,idx]=ismember(txt,fealabels);
        focalg = fealabels(idx);
    end

    function i_genecards(~, ~)
        web(sprintf('https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s', focalg),'-new');
    end
end


