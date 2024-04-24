function [X, genelist, celllist, ftdone, answer1] = i_inputgeolink_mtx

X = [];
genelist = [];
celllist = [];
ftdone = false;
answer1 = [];

prompt = {'Enter link to matrix.mtx.gz:', ...
    'Enter link to genes.tsv.gz:', ...
    'Enter link to barcodes.tsv.gz (OPTIONAL):'};
dlgtitle = 'Input Download Links';
dims = [1, 100];
definput = {'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE150nnn/GSE150741/suppl/GSE150741_HEL24_3_matrix.mtx.gz', ...
    'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE150nnn/GSE150741/suppl/GSE150741_HEL24_3_genes.tsv.gz', ...
    'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE150nnn/GSE150741/suppl/GSE150741_HEL24_3_barcodes.tsv.gz'};
answer = inputdlg(prompt, dlgtitle, dims, definput);

if isempty(answer), return; end
answer1 = strtrim(answer{1});
answer2 = strtrim(answer{2});
answer3 = strtrim(answer{3});

if ~isempty(answer1) && ~isempty(answer2)
    pw1 = pwd();
    a = tempname;
    mkdir(a);
    cd(a);
    fw = gui.gui_waitbar;
    try
        fprintf('Downloading %s...\n', answer2);
        %gunzip(answer{2});
        websave('genes.tsv.gz', answer2);
        fprintf('Downloading %s...\n', answer1);
        %gunzip(answer{1});
        websave('matrix.mtx.gz', answer1);
        if ~isempty(answer3)
            % gunzip(answer{3});
            fprintf('Downloading %s...\n', answer3);
            websave('barcodes.tsv.gz', answer3);
        end
        [X, genelist, celllist, ftdone] = sc_read10xdir(a);
    catch ME
        cd(pw1);
        gui.gui_waitbar(fw, true);
        errordlg(ME.message);
        return;
    end
    gui.gui_waitbar(fw);
    cd(pw1);
end
end


%{
function [X, genelist, celllist, ftdone] = i_inputgeolink_mtx2

FigWidth = 175;
FigHeight = 100;
FigPos(3:4) = [FigWidth, FigHeight]; %#ok

X = [];
genelist = [];
ftdone = [];
celllist = [];

fig = dialog('Position', [300, 300, 540, 150], 'Name', 'My Dialog');

txt1 = uicontrol('Parent', fig, ...
    'Style', 'text', ...
    'Position', [0, 100, 210, 40], ...
    'String', 'Enter link to matrix.mtx.gz:');

txt2 = uicontrol('Parent', fig, ...
    'Style', 'text', ...
    'Position', [0, 50, 210, 40], ...
    'String', 'Enter link to features.tsv.gz:');

btn = uicontrol('Parent', fig, ...
    'Position', [85, 20, 70, 25], ...
    'String', 'Close', ...
    'Callback', 'delete(gcf)');

definput = {'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE150nnn/GSE150741/suppl/GSE150741_HEL24_3_matrix.mtx.gz', ...
    'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE150nnn/GSE150741/suppl/GSE150741_HEL24_3_genes.tsv.gz', ...
    'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE150nnn/GSE150741/suppl/GSE150741_HEL24_3_barcodes.tsv.gz'};

c1 = uicontrol(fig, 'Style', 'edit');
c1.Position = [20, 100, 500, 20];
c1.String = definput(1);
c1.Callback = @change2;

c2 = uicontrol(fig, 'Style', 'edit');
c2.Position = [20, 55, 500, 20];
c2.String = definput(2);
uicontrol(c1);


    function change2(~, ~)
        newc2 = strrep(get(c1, 'String'), 'matrix.mtx', 'genes.tsv');
        set(c2, 'String', newc2);
end

    return;


    prompt = {'Enter link to matrix.mtx.gz:', ...
        'Enter link to genes.tsv.gz:', ...
        'Enter link to barcodes.tsv.gz (OPTIONAL):'};
    dlgtitle = 'Input Download Links';
    dims = [1, 100];

    answer = inputdlg(prompt, dlgtitle, dims, definput);

    if isempty(answer), return; end
    answer1 = strtrim(answer{1});
    answer2 = strtrim(answer{2});
    answer3 = strtrim(answer{3});

    if ~isempty(answer1) && ~isempty(answer2)
        pw1 = pwd();
        a = tempname;
        mkdir(a);
        cd(a);
        fw = gui.gui_waitbar;
        try
            fprintf('Downloading %s...\n', answer2);
            %gunzip(answer{2});
            websave('genes.tsv.gz', answer2);
            fprintf('Downloading %s...\n', answer1);
            %gunzip(answer{1});
            websave('matrix.mtx.gz', answer1);
            if ~isempty(answer3)
                % gunzip(answer{3});
                fprintf('Downloading %s...\n', answer3);
                websave('barcodes.tsv.gz', answer3);
            end
            [X, genelist, celllist, ftdone] = sc_read10xdir(a);
        catch ME
            cd(pw1);
            gui.gui_waitbar(fw, true);
            errordlg(ME.message);
            return;
        end
        gui.gui_waitbar(fw);
        cd(pw1);
    end
end
%}