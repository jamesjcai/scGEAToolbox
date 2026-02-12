function genes = i_str2genelist(s)
%geneListFromString Split a comma-separated gene string into a string array.

    s = string(s);
    s = strip(s);
    s = strip(s, ",");          % remove leading/trailing commas
    genes = split(s, ",");
    genes = strip(genes);
    genes = genes(genes ~= ""); % remove empty entries
end