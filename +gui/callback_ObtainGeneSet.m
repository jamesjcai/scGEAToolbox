function callback_ObtainGeneSet(~,~)

[selecteditem] = gui.i_selectcellscore;
if isempty(selecteditem), return; end
    switch selecteditem
        case 'MSigDB Signature Score...'
            stag=gui.i_selectspecies(2,true);
            if isempty(stag), return; end    
            try
                [posg,ctselected]=gui.i_selectMSigDBGeneSet(stag);
            catch ME
                errordlg(ME.message);
                return;
            end
            if isempty(posg) || isempty(ctselected), return; end
            
            posg
            ttxt = ctselected;

        case 'PanglaoDB Cell Type Marker Score...'
            stag=gui.i_selectspecies(2,true);
            if isempty(stag), return; end
        
            oldpth=pwd;
            pw1=fileparts(mfilename('fullpath'));
            pth=fullfile(pw1,'..','+run','thirdparty','alona_panglaodb');
            cd(pth);
            
            markerfile=sprintf('marker_%s.mat',stag);
            if exist(markerfile,'file')
                load(markerfile,'Tm');
            else
                % Tw=readtable(sprintf('markerweight_%s.txt',stag));
                Tm=readtable(sprintf('markerlist_%s.txt',stag),...
                    'ReadVariableNames',false,'Delimiter','\t');
                % save(markerfile,'Tw','Tm');
            end
            cd(oldpth);
            
            ctlist=string(Tm.Var1);
            listitems=sort(ctlist);
            [indx,tf] = listdlg('PromptString',...
            {'Select Class'},...
             'SelectionMode','single','ListString',listitems,'ListSize',[220,300]);
            if ~tf==1, return; end
            ctselected=listitems(indx);
            % idx=find(matches(ctlist,ctselected));
            idx=matches(ctlist,ctselected);
            ctmarkers=Tm.Var2{idx};
            posg=string(strsplit(ctmarkers,','));
            posg(strlength(posg)==0)=[];
            ttxt = ctselected;
                        
        case 'Select a Predefined Score...'
            gui.gui_showrefinfo('Other Predefined Cell Score');
            [~,T]=pkg.e_cellscores(sce.X,sce.g,0);
            listitems=T.ScoreType;
            [indx2,tf2] = listdlg('PromptString','Select Score',...
                 'SelectionMode','single','ListString',...
                 listitems,'ListSize',[320,300]);
            if tf2~=1, return; end
            [~,~,posg]=pkg.e_cellscores(sce.X,sce.g,indx2);
            ttxt=T.ScoreType(indx2);
        case {'TF Activity Score [PMID:33135076] 🐢',...
                'TF Targets Expression Score...'}
            [~,T]=pkg.e_tfactivityscores(sce.X,sce.g,0);
            listitems=unique(T.tf);

            %[glist]=gui.i_selectngenes(string(listitems));

            [indx2,tf2] = listdlg('PromptString','Select a transcription factor (TF)',...
                 'SelectionMode','single','ListString',...
                 listitems,'ListSize',[220,300]);
            if tf2~=1, return; end

                [cs,tflist]=sc_tfactivity(sce.X,sce.g,[], ...
                    species,methodid);
                idx=find(tflist==string(listitems{indx2}));
                assert(length(idx)==1)
                
                [y]=cs(idx,:);
                ttxt=listitems{indx2};
                posg=[];    % xxxxxxxxxxx

        otherwise
            return;
    end

 % assignin('base','y',y);
 % assignin('base','thisc',thisc);

 if ~exist("posg","var"), posg=[]; end
 if ~exist("ttxt","var"), ttxt=[]; end


end
