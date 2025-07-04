function [TbpUp, TmfUp, TbpDn, TmfDn] = in_XLSX2DETable(excelfile)
    TbpUp = [];
    TmfUp = [];
    TbpDn = [];
    TmfDn = [];
    sheetList = sheetnames(excelfile);
    
    sheetToRead = 'Up_250_GO_BP';
    if any(strcmp(sheetList, sheetToRead))
        TbpUp = readtable(excelfile, 'Sheet', sheetToRead);
    end

    sheetToRead = 'Up_250_GO_MF';
    if any(strcmp(sheetList, sheetToRead))
        TmfUp = readtable(excelfile, 'Sheet', sheetToRead);
    end

    sheetToRead = 'Dn_250_GO_BP';
    if any(strcmp(sheetList, sheetToRead))
        TbpDn = readtable(excelfile, 'Sheet', sheetToRead);
    end

    sheetToRead = 'Dn_250_GO_MF';
    if any(strcmp(sheetList, sheetToRead))
        TmfDn = readtable(excelfile, 'Sheet', sheetToRead);
    end        
end