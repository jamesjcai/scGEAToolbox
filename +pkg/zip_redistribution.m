if ispc
    MyFolderInfo = dir('scGEApp_standalone\for_redistribution_files_only\*');
    fdname = MyFolderInfo(3).folder;
    flname = {''};
    for k = 3:length(MyFolderInfo)
        flname{k-2} = fullfile(fdname, MyFolderInfo(k).name);
    end
    zip('scGEApp_standalone\for_redistribution\scGEAppInstaller.zip', flname);
    web('docs.genomezoo.net', '-browser');
elseif isunix
    MyFolderInfo = dir('scGEApp_standalone/for_redistribution_files_only/*');
    fdname = MyFolderInfo(3).folder;
    flname = {''};
    for k = 3:length(MyFolderInfo)
        flname{k-2} = fullfile(fdname, MyFolderInfo(k).name);
    end
    tar('scGEApp_standalone/for_redistribution/scGEAppInstaller.tar', flname);
    gzip('scGEApp_standalone/for_redistribution/scGEAppInstaller.tar');
    recycle('on');
    delete('scGEApp_standalone/for_redistribution/scGEAppInstaller.tar');
    web('docs.genezoo.org', '-browser');
end
