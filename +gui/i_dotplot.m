%function i_dotplot(X0,X1,genelist,tgene,uselog)
function [hFig]=i_dotplot(X,g,c,cL,tgene,uselog)

if nargin<6, uselog=false; end
[yes]=ismember(tgene,g);
if ~any(yes), return; end
z=length(tgene)-sum(yes);
if z>0
    fprintf('%d gene(s) not in the list are excluded.\n',z); 
end
tgene=tgene(yes);

% tgene=string(T.gene(1:10));
%idx=(1:length(tgene))';
%x=[-ones(size(idx)); ones(size(idx))]./2;
%y=repmat(idx,length(cL),1);

l=ones(length(tgene)*length(cL),1);
sz=l; vl=l;
x=l; y=l;
ct=0;
for kx=1:length(tgene)
    for kk=1:length(cL)
        ct=ct+1;
        x(ct)=kk; y(ct)=kx;
        a0=X(g==tgene(kx),c==kk);
        sz(ct)=sum(a0~=0)./length(a0);
        vl(ct)=mean(a0);
    end
end

% qx=quantile(vl(:),0.10);
% vl(vl<qx)=0;


AvgExpr=vl;
PrtExpr=sz;

GroupList=repmat(string(cL),length(tgene),1);
GeneList=[];
for kx=1:length(tgene)
    GeneList=[GeneList; repmat(tgene(kx),length(cL),1)];
end

% GroupList=[];
% GeneList=repmat(tgene,length(cL),1);
% for k=1:length(cL)
%     GroupList=[GroupList; repmat(string(cL{k}),length(tgene),1)];
% end

% 
% assignin("base","tgene",tgene);
% assignin("base","sz",sz);
% assignin("base","vl",vl);
% assignin("base","cL",cL);

T=table(GeneList,GroupList,AvgExpr,PrtExpr);
% assignin("base","T",T);


sz(sz<0.05)=0;
if uselog
    vl=log2(vl+1);
end
txgene=[" "; tgene];

% figure;
%sz=randi(100,1,length(x));
%scatter([-.5 .5],[-1 -1],[1 500],'k','filled');
%hold on
hFig=figure('Visible','off');

dotsz=1.0;
sz(sz==0)=eps;
vl=vl+0.001;
afa=scatter(x,y,dotsz*500*sz,vl,'filled');
hold on
afb=scatter(x,y,dotsz*500*sz,'k');

af{1}=scatter(max(x)+1,1,dotsz*500*1,'k');
af{2}=text(max(x)+1.4,1,'100%','BackgroundColor','none');
af{3}=scatter(max(x)+1,2,dotsz*500*0.5,'k');
af{4}=text(max(x)+1.4,2,'50%','BackgroundColor','none');
af{5}=scatter(max(x)+1,3,dotsz*500*0.1,'k');
af{6}=text(max(x)+1.4,3,'10%','BackgroundColor','none');


xlim([0.5 length(cL)+2.5]);
ylim([0.5 max([4 length(txgene)])-0.5]);
colorbar
set(gca,'YTick',0:length(tgene))
set(gca,'YTickLabel',txgene)
set(gca,'XTick',0:length(cL))
set(gca,'XTickLabel',[{''};cL(:);{''}])
colormap(flipud(bone));
box on
grid on
% hFig=gcf;
hFig.Position(3)=hFig.Position(3)*0.7;

tb = uitoolbar('Parent', hFig);
pkg.i_addbutton2fig(tb,'on',{@gui.i_pickcolormap,c},'plotpicker-compass.gif','Pick new color map...');
pkg.i_addbutton2fig(tb,'off',@gui.i_changefontsize,'noun_font_size_591141.gif','ChangeFontSize');
pkg.i_addbutton2fig(tb,'off',@i_resizedot,'networkcomp.gif','Resize dots...');
pkg.i_addbutton2fig(tb,'on',@i_renamecat,'guideicon.gif','Rename groups...');
pkg.i_addbutton2fig(tb,'on',{@gui.i_savemainfig,3},"powerpoint.gif",'Save Figure to PowerPoint File...');
pkg.i_addbutton2fig(tb,'on',@i_savetable,'export.gif','Export data...');
pkg.i_addbutton2fig(tb,'on',@gui.i_invertcolor,'plotpicker-comet.gif','Invert colors');
pkg.i_addbutton2fig(tb,'off',@i_resetcolor,'plotpicker-geobubble2.gif','Reset color map');

movegui(hFig, 'center');
set(hFig, 'visible', 'on');

%     function i_changefontsize(~,~)
%           ax=get(gca,'FontSize')+1;
%          if ax>15, ax=5; end
%          set(gca,'FontSize',ax);
%     end

    function i_savetable(~,~)
            answer = questdlg('Export & save data to:','',...
                'Workspace','TXT/CSV file','Excel file','Workspace');
            if ~isempty(answer)
                GroupList=repmat(string(cL),length(tgene),1);
                GeneList=[];
                for k=1:length(tgene)
                    GeneList=[GeneList; repmat(tgene(k),length(cL),1)];
                end
                T=table(GeneList,GroupList,AvgExpr,PrtExpr);
                switch answer
                    case 'Workspace'
                        labels = {'Save T to variable named:'}; 
                        vars = {'T'};
                        values = {T};
                        [~,OKPressed]=export2wsdlg(labels,vars,values,...
                            'Save Data to Workspace');                    
                    case 'TXT/CSV file'
                        [file, path] = uiputfile({'*.csv';'*.*'},'Save as');
                        if isequal(file,0) || isequal(path,0)
                           return;
                        else
                           fw=gui.gui_waitbar; 
                           filename=fullfile(path,file);
                           writetable(T,filename,'FileType','text');
                           OKPressed=true;
                           gui.gui_waitbar(fw);
                        end
                    case 'Excel file'
    
                        [file, path] = uiputfile({'*.xlsx';'*.*'},'Save as');
                        if isequal(file,0) || isequal(path,0)
                           return;
                        else
                           fw=gui.gui_waitbar; 
                           filename=fullfile(path,file);
                           writetable(T,filename,'FileType','spreadsheet');
                           OKPressed=true;
                           gui.gui_waitbar(fw);
                        end
                end
            end
    end

    function i_resizedot(~,~)
        dotsz=dotsz*0.9;
        if dotsz<0.2, dotsz=1.0; end
        delete(afa); delete(afb);
        delete(af{1}); delete(af{3}); delete(af{5});
        afa=scatter(x,y,dotsz*500*sz,vl,'filled');
        hold on
        afb=scatter(x,y,dotsz*500*sz,'k');
        af{1}=scatter(max(x)+1,1,dotsz*500*1,'k');
        af{3}=scatter(max(x)+1,2,dotsz*500*0.5,'k');
        af{5}=scatter(max(x)+1,3,dotsz*500*0.1,'k');
    end

    function i_renamecat(~,~)
        tg=gui.i_inputgenelist(string(cL),true);
        if isempty(tg), return; end
        if length(tg)==length(cL)
            set(gca,'XTick',0:length(cL));
            set(gca,'XTickLabel',[{''};tg(:);{''}])
            cL=tg;
        else
            errordlg('Wrong input.');
        end        
    end

    function i_resetcolor(~,~)
        dotsz=0.1;
        set(gca,'FontSize',10);
        i_resizedot;
        colormap(flipud(bone));
    end

%     function i_invertcolor(~,~)
%         cm=colormap();
%         colormap(flipud(cm));
%     end
end
