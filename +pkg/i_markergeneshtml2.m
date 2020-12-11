function i_markergeneshtml2(sce,markerlist,numfig,axbx)
if nargin<4, axbx=[]; end
if nargin<3, numfig=20; end

numfig=min([numfig length(markerlist)]);
dirtxt=tempdir;
a=colormap('autumn');
a(1,:)=[.8 .8 .8];

fname=[tempname,'.html'];
fid=fopen(fname,'w');
htmlstr="";
    for k=1:numfig
        k
        targeetg=markerlist(k);
        h=figure('Visible','off');
%         sc_scattermarker(sce.X,sce.g,sce.s,...
%             markerlist(k),3,5,false);        
        z=log2(1+sce.X(sce.g==targeetg,:));
        subplot(1,2,1)
                scatter(sce.s(:,1),sce.s(:,2),...
                1,z,'filled');
            box on
%         scatter3(sce.s(:,1),sce.s(:,2),sce.s(:,3),...
%                 5,z,'filled');
        colormap(h,a);
        title(targeetg)
        if ~isempty(axbx), view(axbx(1),axbx(2)); end
        %imgfname1=sprintf('heatmap_%s.png',targeetg);
        %saveas(h,sprintf('%s%s',dirtxt,imgfname1));
        %close(h);
        
        %h=figure('Visible','off');
        subplot(1,2,2)
        
[~,cL]=grp2idx(sce.c);
[~,i]=sort(grpstats(z,sce.c,@median),'descend');
xn=matlab.lang.makeValidName(cL{i(1)});
 
        pkg.i_violinplot_groupordered(z,sce.c);
        ylabel('log2(1+UMI)');
        title(targeetg);
        xtickangle(-45);
        %suptitle(h,targeetg);
        
        imgfname=sprintf('%s.png',targeetg);
        
        targetdir=sprintf('G:\\My Drive\\Ahr differentiation potency\\marker_candidates\\%s',xn);
        if ~exist(targetdir,'dir')
            mkdir(targetdir);
        end
        saveas(h,sprintf('%s%s',targetdir,imgfname));
        close(h);
        
%         aax=sprintf('<table>\n<tr>\n<th><img src="%s"></th>\n</tr>\n</table>\n',...
%             imgfname);
        aax=sprintf('<img src="%s"><br>\n',imgfname);
        %htmlstr=sprintf('%s\n%s',htmlstr,aax);
    end
%    fprintf(fid,'%s',htmlstr);
    fclose(fid);
    pause(1);
    % winopen(tempdir);
    web(fname,'-browser');
end
