function callback_scTenifoldXct(src,~)
    % import ten.*
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);

if ~gui.i_setpyenv
    return; 
end

    [thisc,~]=gui.i_select1class(sce,false);
    if isempty(thisc), return; end    
    [c,cL]=grp2idx(thisc);
    [idx]=gui.i_selmultidlg(cL);
    if isempty(idx), return; end
    if numel(idx)<2
        warndlg('Need at least 2 cell groups to perform cell-cell interaction analysis.');
        return;
    end
    if numel(idx)~=2
        warndlg(sprintf('Need only 2 cell groups to perform cell-cell interaction analysis. You selected %d.',...
                numel(idx)));
        return;
    end
    
    i1=idx(1);
    i2=idx(2);


%{
    
[~,cL]=grp2idx(sce.c_cell_type_tx);
if length(cL)<2, errordlg('Need at least 2 cell types.'); return; end

[indxx,tf2] = listdlg('PromptString',...
    {'Select two cell types:'},...
     'SelectionMode','multiple','ListString',cL);
if tf2==1
    if numel(indxx)~=2
        errordlg('Please select 2 cell types');
        return;
    end
    i1=indxx(1);
    i2=indxx(2);
else
    return;
end
%}


a1=sprintf('%s -> %s',cL{i1},cL{i2});
a2=sprintf('%s -> %s',cL{i2},cL{i1});

twosided=false;
answer=questdlg('Select direction: Source (ligand) -> Target (receptor)','','Both',a1,a2,'Both');
switch answer
    case 'Both'
        x1=i1; x2=i2;
        twosided=true;
    case a1
        x1=i1; x2=i2;        
    case a2        
        x1=i2; x2=i1;        
    otherwise
        return;
end

%{
idx=sce.c_cell_type_tx==cL{x1} | sce.c_cell_type_tx==cL{x2};
sce=sce.selectcells(idx);

sce.c_batch_id=sce.c_cell_type_tx;
sce.c_batch_id(sce.c_cell_type_tx==cL{x1})="Source";
sce.c_batch_id(sce.c_cell_type_tx==cL{x2})="Target";
%}



sce.c_batch_id=thisc;
sce.c_batch_id(c==x1)="Source";
sce.c_batch_id(c==x2)="Target";
sce.c_cell_type_tx=string(cL(c));

% idx=thisc==cL{x1} | thisc==cL{x2};
idx=c==x1 | c==x2;
sce=sce.selectcells(idx);


%sce.c_batch_id(thisc==cL{x1})="Source";
%sce.c_batch_id(thisc==cL{x2})="Target";

try    
    if twosided
        [Tcell]=run.py_scTenifoldXct(sce,cL{x1},cL{x2},true);
        [T1]=Tcell{1};
        [T2]=Tcell{2};
        if ~isempty(T1)
            a=sprintf('%s -> %s',cL{x1},cL{x2});
            T1 = addvars(T1,repelem(a,size(T1,1),1),'Before',1);
            T1.Properties.VariableNames{'Var1'} = 'direction';
        end
        if ~isempty(T2)      
            a=sprintf('%s -> %s',cL{x2},cL{x1});
            T2 = addvars(T2,repelem(a,size(T2,1),1),'Before',1);
            T2.Properties.VariableNames{'Var1'} = 'direction';
        end
        T=[T1;T2];
    else
        [T]=run.py_scTenifoldXct(sce,cL{x1},cL{x2},false);
        %T=readtable('output1.txt');
        if ~isempty(T)
            a=sprintf('%s -> %s',cL{x1},cL{x2});
            T = addvars(T,repelem(a,size(T,1),1),'Before',1);
            T.Properties.VariableNames{'Var1'} = 'direction';            
        end
    end
catch ME
    errordlg(ME.message);
    return;
end

if ~isempty(T)
    [b,a]=pkg.i_tempfile("sctendifoldxct");
    writetable(T,b);

    answer=questdlg(sprintf('Result has been saved in %s',b), ...
        '','Export result...','Locate result file...','Export result...');
    switch answer
        case 'Locate result file...'
            winopen(a);
            pause(2)
            reshowdlg;
        case 'Export result...'
            gui.i_exporttable(T);
        otherwise
            winopen(a);            
    end    
else
    helpdlg('No ligand-receptor pairs are identified.','');
end


    function reshowdlg
        answer=questdlg('Export result to other format?','');
        switch answer
            case 'Yes'
                gui.i_exporttable(T);
            otherwise
                return;
        end        
    end

end