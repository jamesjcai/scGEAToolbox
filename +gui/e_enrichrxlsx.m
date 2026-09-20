function e_enrichrxlsx(Tup,Tdn,T,filesaved)
% Sheets are written only for libraries that return terms. A list the DE
% filter left empty, or one too short to enrich anything, used to leave
% the workbook without enrichment sheets and without saying why - which
% reads as a failure of the enrichment step rather than of the filter.
in_enrichdirection(Tup, T, filesaved, 'Up');
in_enrichdirection(Tdn, T, filesaved, 'Dn');
end

    function in_enrichdirection(Tdir, T, filesaved, dirtag)
        libraries = ["GO_Biological_Process_2025", ...
                     "GO_Molecular_Function_2025", ...
                     "KEGG_2021_Human", ...
                     "Reactome_Pathways_2024"];
        sheettags = ["GO_BP", "GO_MF", "KEGG", "Reactome"];
        [~, fname, fext] = fileparts(filesaved);
        fname = string(fname) + string(fext);

        ngene = height(Tdir);
        if ngene == 0
            fprintf(['%s: %s list has 0 genes after DE filtering - ' ...
                'enrichment skipped. Relax the DE filter to enrich it.\n'], ...
                fname, dirtag);
            return;
        end

        Tlist = run.ml_Enrichr(Tdir.gene(1:min([250 ngene])), T.gene, libraries);
        nterm = cellfun(@height, Tlist);
        for k = 1:numel(libraries)
            in_writetable(Tlist{k}, filesaved, ...
                sprintf('%s_250_%s', dirtag, sheettags(k)));
        end
        if all(nterm == 0)
            fprintf(['%s: %s list (n = %d) returned no enriched terms - ' ...
                'no enrichment sheets written.\n'], fname, dirtag, ngene);
        end
    end

    function in_writetable(Tmf1, filesaved, shtname)
        if ~isempty(Tmf1) && istable(Tmf1) && height(Tmf1) > 0
            if isExcelFile(filesaved)
                writetable(Tmf1, filesaved, "FileType", "spreadsheet", ...
                    'Sheet', shtname);
            else
                filename = "enrichrxlsx_" + matlab.lang.makeValidName(shtname) + ".txt";
                writetable(Tmf1, filename, "FileType", "text");
            end
        end
    end


    function isExcel = isExcelFile(filename)
        isExcel = endsWith(filename, {'.xls', '.xlsx'}, 'IgnoreCase', true);
        %{
        isExcel = contains(filename, {'.xls', '.xlsx'}, 'IgnoreCase', true);
        [~, ~, ext] = fileparts(filename);
        if ~ismember(lower(ext), {'.xls', '.xlsx'})
            isExcel = false;
            return;
        end

        fid = fopen(filename, 'r');
        if fid == -1
            isExcel = false;
            return;
        end
        bytes = fread(fid, 4, 'uint8')';
        fclose(fid);

        % Excel file signatures
        isExcel = isequal(bytes, [208 207 17 224]) || ... % .xls (OLE Compound File)
                  isequal(bytes, [80 75 3 4]);          % .xlsx (ZIP-based)
        %}
    end
