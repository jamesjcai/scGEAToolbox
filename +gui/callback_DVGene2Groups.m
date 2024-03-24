function callback_DVGene2Groups(src, ~)

    FigureHandle = src.Parent.Parent;
    sce = guidata(FigureHandle);
    gui.gui_showrefinfo('DV Analysis');
    
    
    [i1, i2, cL1, cL2] = gui.i_select2grps(sce);
    if length(i1) == 1 || length(i2) == 1, return; end
    
    c=zeros(size(i1));
    c(i1)=1; c(i2)=2;
    cL=[cL1;cL2];
    if ~all(c>0)
        sce=sce.selectcells(c>0);
        c=c(c>0);
        i1=c==1;
        i2=c==2;
    end
    
    % fw = gui.gui_waitbar;
    
    [T1, ~, ~, xyz1] = sc_splinefit(sce.X(:,i1), sce.g, false, false, true);
    [T2, ~, ~, xyz2] = sc_splinefit(sce.X(:,i2), sce.g, false, false, true);
    
    valididx = T1.nearidx>0 & T2.nearidx>0;
    
    bd = abs(xyz1(T1.nearidx(valididx)) - xyz2(T2.nearidx(valididx))); % baseline difference
    dd = abs(T1.d(valididx) - T2.d(valididx));
    ddn = dd./bd;
    
    BaselineDiffDist = zeros(size(T1,1), 1);
    DiffDistRaw =  zeros(size(T1,1), 1);
    DiffDistNormlized = zeros(size(T1,1), 1);
    
    BaselineDiffDist(valididx) = bd;
    DiffDistRaw(valididx) = dd;
    DiffDistNormlized(valididx) = ddn;
    
    T1.Properties.VariableNames = append(T1.Properties.VariableNames, sprintf('_%s',cL1{1}));
    T2.Properties.VariableNames = append(T2.Properties.VariableNames, sprintf('_%s',cL2{1}));
    
    T=[T1 T2 table(BaselineDiffDist, DiffDistRaw, DiffDistNormlized)];
    
    %assignin("base","T",T)
    %assignin("base","T1",T1)
    %assignin("base","T2",T2)
    %assignin("base","xyz1",xyz1)
    %assignin("base","xyz2",xyz2)
    
    T = sortrows(T,"DiffDistNormlized","descend");
    
    %    gui.gui_waitbar(fw);

    outfile = sprintf('%s_vs_%s', ...
        matlab.lang.makeValidName(string(cL1)), matlab.lang.makeValidName(string(cL2)));
    [~, filesaved] = gui.i_exporttable(T, true, 'Tdvgenelist', outfile);
    if ~isempty(filesaved)
       waitfor(helpdlg(sprintf('Result has been saved in %s',filesaved),''));
    end
end