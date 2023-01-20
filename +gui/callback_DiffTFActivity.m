function callback_DiffTFActivity(src,~)
    
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);


    [thisc,clable]=i_select1class(sce,false);

%     if length(unique(thisc))>2
%         answer = questdlg('Two group or multiple group?');
%         switch answer
%             case 'Yes'
%                 [i1,i2,cL1,cL2]=i_select2grps(thisc);
%                 twogroup=true;
%             case 'No'      
%                 twogroup=false;
%             otherwise
%                 return;
%         end
%     end

twogroup=false;

    species=gui.i_selectspecies;

    [cs,tflist,gcommon]=sc_tfactivity(sce.X,sce.g,[],species);
    idx=grp2idx(thisc);
    
    P=ones(length(tflist),1);
    D=true(length(tflist),1);
    for k=1:length(tflist)
        s=cs(k,:);
        if twogroup
            %[~,p]=kstest2(s(idx==1),s(idx==2));
            [~,p]=ttest2(s(idx==1),s(idx==2));
            if median(s(idx==1))>median(s(idx==2))
                D(k)=false;
            end            
        else
            p = anovan(s,idx);            
        end
        P(k)=p;
    end
    T=table(tflist,P);
    T=sortrows(T,"P","ascend");

%     outfile=sprintf('%s_vs_%s', ...
%         matlab.lang.makeValidName(string(cL1)),matlab.lang.makeValidName(string(cL2)));    
      outfile=sprintf('DiffTFActivity_%s',matlab.lang.makeValidName(clable));

    [filetype,filesaved]=gui.i_exporttable(T,true,'T',outfile);

end