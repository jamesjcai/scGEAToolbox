function [X,genelist,celllist,ftdone]=i_inputgeolinks

    X=[]; genelist=[]; ftdone=[]; celllist=[];
    prompt = {'Enter link to matrix.mtx.gz:',...
        'Enter link to genes.tsv.gz:',...
        'Enter link to barcodes.tsv.gz (OPTIONAL):'};
    dlgtitle = 'Input Download Links';
    dims = [1 100];
    definput = {'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE150nnn/GSE150741/suppl/GSE150741_HEL24_3_matrix.mtx.gz',...
        'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE150nnn/GSE150741/suppl/GSE150741_HEL24_3_genes.tsv.gz',...
        'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE150nnn/GSE150741/suppl/GSE150741_HEL24_3_barcodes.tsv.gz'};
    answer = inputdlg(prompt,dlgtitle,dims,definput);
    if isempty(answer), return; end
    if ~isempty(answer{1}) && ~isempty(answer{2})
        pw1=pwd();
        a=tempname;
        mkdir(a);
        cd(a);
        fw=gui.gui_waitbar;
        try
        fprintf('Downloading %s...\n',answer{2});
        %gunzip(answer{2});
        websave('genes.tsv.gz',answer{2});
        fprintf('Downloading %s...\n',answer{1});
        %gunzip(answer{1});
        websave('matrix.mtx.gz',answer{1});
        if ~isempty(answer{3})
            % gunzip(answer{3});
            fprintf('Downloading %s...\n',answer{3});
            websave('barcodes.tsv.gz',answer{3});
        end
        [X,genelist,celllist,ftdone]=sc_read10xdir(a);
        catch ME
            cd(pw1);
            gui.gui_waitbar(fw);
            errordlg(ME.message);
            return;
        end
        gui.gui_waitbar(fw);
        cd(pw1);
    end
end