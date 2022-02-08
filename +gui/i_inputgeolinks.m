function [X,genelist,celllist,ftdone,answer1]=i_inputgeolinks

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
    answer1=strtrim(answer{1});
    answer2=strtrim(answer{2});
    answer3=strtrim(answer{3});    
    
    if ~isempty(answer1) && ~isempty(answer2)
        pw1=pwd();
        a=tempname;
        mkdir(a);
        cd(a);
        fw=gui.gui_waitbar;
        try
        fprintf('Downloading %s...\n',answer2);
        %gunzip(answer{2});
        websave('genes.tsv.gz',answer2);
        fprintf('Downloading %s...\n',answer1);
        %gunzip(answer{1});
        websave('matrix.mtx.gz',answer1);
        if ~isempty(answer3)
            % gunzip(answer{3});
            fprintf('Downloading %s...\n',answer3);
            websave('barcodes.tsv.gz',answer3);
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