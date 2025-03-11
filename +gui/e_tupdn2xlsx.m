function e_tupdn2xlsx(Tup,Tdn,T,filesaved)

    writetable(T, filesaved, 'FileType', 'spreadsheet', 'Sheet', 'All_genes');
    writetable(Tup, filesaved, "FileType", "spreadsheet", 'Sheet', 'Up-regulated');
    writetable(Tdn, filesaved, "FileType", "spreadsheet", 'Sheet', 'Down-regulated');