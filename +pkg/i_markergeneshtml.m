function i_markergeneshtml(sce,markerlist,numfig,axbx,nametag)
if nargin<5, nametag=''; end
if nargin<4, axbx=[]; end
if nargin<3, numfig=20; end

numfig=min([numfig length(markerlist)]);
dirtxt=tempdir;
a=colormap('autumn');
a(1,:)=[.8 .8 .8];

if isempty(nametag)
    htmlfilenamex=tempname;
else
    nametag=matlab.lang.makeValidName(nametag);
    htmlfilenamex=[dirtxt,nametag];
end

fname=[htmlfilenamex,'.html'];
fid=fopen(fname,'w');
htmlstr="";
    for k=1:numfig
        targeetg=markerlist(k);
        h=figure('Visible','off');
%         sc_scattermarker(sce.X,sce.g,sce.s,...
%             markerlist(k),3,5,false);        
        c=log2(1+sce.X(sce.g==targeetg,:));
        scatter3(sce.s(:,1),sce.s(:,2),sce.s(:,3),...
                5,c,'filled');
        colormap(h,a);
        title(targeetg)
        if ~isempty(axbx)
            view(axbx(1),axbx(2));
        end
        imgfname1=sprintf('heatmap_%s.png',targeetg);
        saveas(h,sprintf('%s%s',dirtxt,imgfname1));
        close(h);
        
        h=figure('Visible','off');
        pkg.i_violinplot_groupordered(c,sce.c);
        ylabel('log2(1+UMI)');
        title(targeetg)
        xtickangle(-45);
        imgfname2=sprintf('violin_%s.png',targeetg);
        saveas(h,sprintf('%s%s',dirtxt,imgfname2));
        close(h);
        
        aax=sprintf('<center><table>\n<tr>\n<th><img src="%s"></th>\n<th><img src="%s"></th>\n</tr>\n</table></center>\n',...
            imgfname1,imgfname2);
        htmlstr=sprintf('%s\n%s',htmlstr,aax);
    end
    fprintf(fid,'%s',htmlstr);
    fclose(fid);
    pause(1);
    % winopen(dirtxt);
    web(fname,'-browser');
end
