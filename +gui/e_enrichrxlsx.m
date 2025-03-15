function e_enrichrxlsx(Tup,Tdn,T,filesaved)
    [Tlist1] = run.ml_Enrichr(Tup.gene(1:min([250 height(Tup)])), ...
                T.gene, ["GO_Biological_Process_2025", ...
                         "GO_Molecular_Function_2025", ...
                         "KEGG_2021_Human",...
                         "Reactome_Pathways_2024"]);
    Tbp1 = Tlist1{1};
    Tmf1 = Tlist1{2};
    Keg1 = Tlist1{3};
    Rea1 = Tlist1{4};
    in_writetable(Tbp1, filesaved, 'Up_250_GO_BP');
    in_writetable(Tmf1, filesaved, 'Up_250_GO_MF');
    in_writetable(Keg1, filesaved, 'Up_250_KEGG');
    in_writetable(Rea1, filesaved, 'Up_250_Reactome');

    [Tlist2] = run.ml_Enrichr(Tdn.gene(1:min([250 height(Tdn)])), ...
                T.gene, ["GO_Biological_Process_2025", ...
                         "GO_Molecular_Function_2025",...
                         "KEGG_2021_Human",...
                         "Reactome_Pathways_2024"]);
    Tbp2 = Tlist2{1};
    Tmf2 = Tlist2{2};
    Keg2 = Tlist2{3};
    Rea2 = Tlist2{4};
    in_writetable(Tbp2, filesaved, 'Dn_250_GO_BP');
    in_writetable(Tmf2, filesaved, 'Dn_250_GO_MF');
    in_writetable(Keg2, filesaved, 'Dn_250_KEGG');
    in_writetable(Rea2, filesaved, 'Dn_250_Reactome');
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


