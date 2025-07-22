classdef WebDownload
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu
%   Math Lead & Secondary Developer: Connor Meehan <connor.gw.meehan@gmail.com>
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause

    methods(Static)
        
        function yes=GoogleDriveOnly(yes)
            %helpful for white box testing of Google Drive fail over logic
        
            if nargin>0
                assert(islogical(yes));
                priorYes=WebDownload.Globals('googleDriveOnly', yes);
                app=BasicMap.Global;
                app.urlMap.clear;
                app.set(WebDownload.PROP_GOOGLE_DRIVE_ONLY, num2str(yes));
                yes=priorYes;
            else
                yes=WebDownload.Globals('googleDriveOnly');
            end
            if isempty(yes)
                yes=false;
            end
        end
        
        function priorUrl=UseMeehanGoogleDrive
            priorUrl=WebDownload.GoogleDriveDirectory([]);
        end
        
        function url=GoogleDriveDirectory(url)
            %helpful for white box testing of Google Drive logic
        
            if nargin>0
                if isempty(url)
                    url=WebDownload.MEEHAN_GOOGLE_DRIVE_DIRECTORY;
                end
                assert(ischar(url));
                priorUrl=WebDownload.Globals('googleDriveDirectory', url);
                J=edu.stanford.facs.swing.WebDownload;
                J.SetDefaultGoogleDirUrl(url);
                newFile=fullfile(File.Home, 'myGoogleDir.properties');
                if exist(newFile, 'file')
                    delete(newFile);
                end
                problems=java.util.ArrayList;
                try
                    J.SetDefaultGoogleDirFile(newFile, problems);
                catch ex
                    problems.add(ex.message)
                end
                if problems.size>0
                    try
                        warning('%d problems reading %s\n',...
                            url);
                        it=problems.iterator;
                        while it.hasNext
                            fprintf('%s\n', char(it.next));
                        end
                    catch
                    end
                end
                app=BasicMap.Global;
                app.urlMap.clear;
                app.set(WebDownload.PROP_GOOGLE_DRIVE_DIRECTORY, url);
                url=priorUrl;
            else
                url=WebDownload.Globals('googleDriveDirectory');
                if isempty(url)
                    url=WebDownload.MEEHAN_GOOGLE_DRIVE_DIRECTORY;
                end
            end
        end
        
        function value=Globals(name, value)
            persistent stuff;
            if isempty(stuff)
                stuff=struct();
                app=BasicMap.Global;
                WebDownload.GoogleDriveDirectory(...
                    app.get(...
                    WebDownload.PROP_GOOGLE_DRIVE_DIRECTORY));
                WebDownload.GoogleDriveOnly(...
                    app.is(...
                    WebDownload.PROP_GOOGLE_DRIVE_ONLY, false));
            end
            if nargin>1
                prior=WebDownload.Globals(name);
                stuff.(name) = value;
                value=prior;
            else
                if isfield(stuff, name)
                    value=stuff.(name);
                else
                    value=[];
                end
            end
        end

        function ok=AskDownload(localFile, uri)
            app=BasicMap.Global;
            f1=Html.FileTree(localFile);
            f2=Html.FileTree(uri);
            q=['<html><font color="red">' app.h3Start 'You will ' ...
                'lose your previous version...' app.h3End ...
                '</font><table><thead><tr>' ...
                '<th><b>Overwriting</b></th>' ...
                '<th><b>With download from</b></th></tr>' ...
                '</thead><tr><td>' f1 '</td>' ...
                '<td>' f2 '</td></tr></table>' ...
                app.h2Start '<font color="red">Continue?</font>' ...
                app.h2End '</html>'];
            ok=askYesOrNo(q, 'WARNING:  file overwriting');
        end
    end

    properties(Constant)
        
        %Call SetGoogleDir if the fall Google Drive is not Stephen Meehan's
        MEEHAN_GOOGLE_DRIVE_DIRECTORY='https://drive.google.com/file/d/1r_ByRnrq7SWvfLa2uGdwE31IY1ClXVdy/view?usp=sharing';
        
        PROP_GOOGLE_DRIVE_DIRECTORY='WebDownload.GoogleDriveDirectory';
        PROP_GOOGLE_DRIVE_ONLY='WebDownload.GoogleDriveOnly';
        PATH='GetDown2/demo';
        DOCUMENTS_FOLDER_FOR_EXAMPLES='run_umap/examples';
        
        KEY_UMAP_EXAMPLES='umapExamples';
        KEY_DEMO='demo';
        KEY_AutoGate='facs';
        PROPERTY_IS_FIRST_GOOGLE_CLOUD='IsFirstGoogleCloudDownload';
        PROPERTY_HOST='host';
        DEFAULT_HOST='https://storage.googleapis.com/cytogenie.org';
        HOST=[WebDownload.DEFAULT_HOST '/'];
        
        DEFAULT_BAD_HOST='http://cgworkspace.cytogenies.org';
        HOST_FILE='.herzenbergLabHosts.properties';        
    end
    
    methods(Static)
        
        %GetHosts looks up a list of URL root "mirrors" expected to contain
        %the folder of urlPath.  The internal list is used if this file is
        %not found <HOME>/.herzenbergLabHosts.properties If on the other
        %hand the properties file is found then the first server looked up
        %is indicated by the property host.1 and if this server is
        %unavailable the URL root in host.2 is used. The next server at
        %host.3 is checked if host.2 is also unavailable and the lookup
        %stops when # of hosts specified by host.count is checked.
        %Add or alter the host list as you need to and be sure host.count
        %is correct.
        function [hosts, googleDir, m]=GetHosts(urlPath)
            googleDriveOnly=WebDownload.GoogleDriveOnly;
            if ~googleDriveOnly
                fileName=fullfile(File.Home, WebDownload.HOST_FILE);
                m=JavaProperties(fileName);
                m1=m.addIfMissing(WebDownload.PROPERTY_HOST, WebDownload.DEFAULT_HOST);
                m2=m.addIfMissing(WebDownload.PROPERTY_HOST, 'https://drive.google.com');
                m3=m.addIfMissing(WebDownload.PROPERTY_HOST, 'http://54.39.2.45');
                m4=m.putIfMissing(WebDownload.PROP_GOOGLE_DRIVE_DIRECTORY, ...
                    WebDownload.GoogleDriveDirectory);
                if m1 || m2 || m3 || m4
                    m.save(fileName);
                end
            else
                m=JavaProperties;
                m.set(WebDownload.PROPERTY_HOST, 'https://drive.google.com');
                m.set(WebDownload.PROP_GOOGLE_DRIVE_DIRECTORY, ...
                    WebDownload.GoogleDriveDirectory);
            end
            specificKey=[WebDownload.PROPERTY_HOST '.' urlPath];
            hosts0=m.getAll(specificKey);
            specificKey=[WebDownload.PROP_GOOGLE_DRIVE_DIRECTORY '.' urlPath];
            googleDir=m.get(specificKey);
            hosts=[hosts0 m.getAll(WebDownload.PROPERTY_HOST)];
            if ~googleDriveOnly  && ~isequal(hosts{1}, WebDownload.DEFAULT_HOST)
                app=BasicMap.Global;
                if ~app.is(WebDownload.PROPERTY_IS_FIRST_GOOGLE_CLOUD)
                    %prioritize Google Cloud in local properties if first time
                    hosts=StringArray.Sort(hosts, {WebDownload.DEFAULT_HOST});
                    app.set(WebDownload.PROPERTY_IS_FIRST_GOOGLE_CLOUD, 'yes');
                    app.save;
                    m.removeAll(WebDownload.PROPERTY_HOST);
                    m.setAll(WebDownload.PROPERTY_HOST, hosts)
                    m.save(fileName);
                end
            end
            if isempty(googleDir)
                googleDir=m.get(WebDownload.PROP_GOOGLE_DRIVE_DIRECTORY);
            end
        end
        
        function StateSharingRequirements
            disp('Any URL like https://drive.google.com/run_umap/examples/sample10k.csv');
            disp('  gets converted to an internal Google Drive export format  ');
            disp('  following the logic at https://www.wonderplugin.com/online-tools/google-drive-direct-link-generator/');
            disp('  This conversion requires a properties file ');
            disp('  that records the relative filepath and size and the sharing');
            disp('  link produced by the Google Drive user interface');
            disp('  ');
            disp('THEN code like WebDownload.java can directly download the ');
            disp('  URL in export format ONLY if both the file ');
            disp('  and the parent folder are shared read-only to everyone!!');
            disp(' ');
            disp(' This does not work for files like tar that are too big to download');
        end
        
        function [url, googleDir]=ResolveHost(urlPath)
            app=BasicMap.Global;
            map=app.urlMap;
            [hosts, googleDir]=WebDownload.GetHosts(urlPath);
            N=length(hosts);
            for i=1:N
                beforePath=hosts{i};
                [host,~,port]=WebDownload.UrlParts(beforePath);
                if isempty(port) || port<=0 
                    port=80;
                end
                badKey=['bad:' host ':' num2str(port)];
                if ~map.containsKey(badKey)
                    [ok, ~]=WebDownload.CanReach(host, port);
                    if ok
                        url=WebDownload.FullUrl(beforePath, urlPath);
                        return;
                    end
                    map.put(badKey, badKey);
                end
            end
            if N>0
                map.clear;
                url=urlPath;
            else
                url='';
            end
        end
        
        function [ok, issues]=CanReach(host, port, timeout)
            if nargin<3
                timeout=4000;
            end
            if isempty(port) || port<=0 
                port=80;
            end
            issues=java.util.ArrayList;
            ok=edu.stanford.facs.swing.WebDownload.CanReachHost(host, port, timeout, issues);
        end
        
        function [url, googleDir]=ResolveUrl(file, key, specialHost)
            if nargin<2
                key='run_umap/examples';
                if nargin<1
                    file='';
                end
            end
            app=BasicMap.Global;
            map=app.urlMap;
            firstTime=map.size==0;
            googleDir=[];
            if WebDownload.GoogleDriveOnly
                map.clear;
            end
            if nargin>2 && ~isempty(specialHost)
                if endsWith(specialHost, '/')
                    url=[specialHost key];
                else
                    url=[specialHost '/' key];
                end
            elseif ~map.containsKey(key)
                [url, googleDir]=WebDownload.ResolveHost(key);
                map.put(key, url);
            else
                url=map.get(key);
            end
            if ~firstTime
                [host,~,port]=WebDownload.UrlParts(url);
                if ~WebDownload.CanReach(host, port)
                    map.clear;
                    [url, googleDir]=WebDownload.ResolveUrl(file, key);
                end
            end
            if ~isempty(url)
                url=[url '/' strrep(file, '\', '/')];
            end
        end
        
        function url=Url(file, path, host)
            if nargin<3
                host=WebDownload.HOST;
                if nargin<2
                    path=WebDownload.PATH;
                    if nargin<1
                        file='';
                    end
                end
            end
            url=[host path '/' file];
        end
        
        function [fullFile, exists, downloaded]=Download(...
                fromUrl, whereOnScreen, toLocalFolder, allowCancel)
            if nargin<4
                allowCancel=true;
                if nargin<3
                    toLocalFolder=tempdir;
                    if nargin<2
                        whereOnScreen='center';         
                    end
                end
            end
            exists=false;
            downloaded=false;
            [host,path]=WebDownload.UrlParts(fromUrl);
            if isempty(host)
                fullFile=[];
                msgError(Html.WrapHr(['Not a URL...<br>'...
                    Html.WrapSmallTags(['(Bad format "' ...
                    fromUrl '")' ])]));
                return;
            end
            if isempty(path)
                remoteFile='tempWebDownloadFile';
            else
                [~, f, e]=fileparts(path);
                remoteFile=[f e];
            end
            fullFile=fullfile(toLocalFolder, remoteFile);
            if exist(fullFile, 'file') && ~isequal(tempdir, toLocalFolder)
                if ~askYesOrNo(Html.WrapHr(['Overwrite prior ?<br>' ...
                        Html.WrapSmallTags(['(' fullFile ')'])]), ...
                        'Exists  ....', whereOnScreen, true)
                    exists=true;
                    return;
                end
            end
            [cancelledByUser, bad]=WebDownload.Get(...
                {fromUrl}, {fullFile}, false, allowCancel, whereOnScreen); 
            downloaded=~cancelledByUser && ~bad;
            fprintf('Downloaded=%d, cancelled=%d and bad=%d\n', ...
                downloaded, cancelledByUser, bad);
            if ~downloaded
                fullFile=[];
            end                
        end
        
        function [bytesPerUrl, didUrlConnect, problemsPerUrl, dwl]=GetSize(urls)
            try
                dwl=javaObjectEDT('edu.stanford.facs.swing.WebDownload');
            catch
                jar=fullfile(fileparts(mfilename('fullpath')), 'webDownload.jar');
                javaaddpath(jar);
                try
                    dwl=javaObjectEDT('edu.stanford.facs.swing.WebDownload');
                catch ex
                    ex.getReport
                    error('Trouble using WebDownload');
                end
            end
            dwl.waitWhenDone=false;
            dwl.allowCancel=false;
            N=length(urls);
            bytesPerUrl=zeros(1,N);
            didUrlConnect=false(1,N);
            problemsPerUrl=cell(1,N);
            for i=1:N
                problemsPerUrl{i}=java.util.ArrayList;
            	con=dwl.connect(urls{i}, problemsPerUrl{i});
                if ~isempty(con.http)
                    didUrlConnect(i)=true;
                    bytesPerUrl(i)= con.size;
                    con.disconnect();
                end		
            end
        end
        
        function [ok, cancelled, dwl, dlg, notExistsFiles,...
                existsFiles, notExistsUrls, existsUrls]=...
                FromUrlFileIfNotExists(urlFile, toFolder)
            try
                lines=strsplit(File.ReadTextFile(urlFile));
                N=length(lines);
                urls=cell(1,N);
                nUrls=0;
                for i=1:N
                    line=lines{i};
                    if isempty(strtrim(line))
                        continue;
                    end
                    [~,f,e]=fileparts(line);
                    file=fullfile(toFolder, String.URLDecode([f e]));
                    if ~exist(file, 'file')
                        nUrls=nUrls+1;
                        urls{nUrls}=line;
                    end
                end
                urls=urls(1:nUrls);
            catch ex
                ex.getReport
                ok=false;cancelled=true;dwl=[];dlg=[];
                notExistsFiles={}; existsFiles={};
                notExistsUrls={};existsUrls={};
                return;
            end
            if nUrls>0
                Gui.Modem;
            end
            [ok, cancelled, dwl, dlg, notExistsFiles, ...
                existsFiles, notExistsUrls, existsUrls]=...
                WebDownload.Many(urls, toFolder);
        end
        
        function [ok, cancelled, dwl, dlg, notExistsFiles,...
                 existsFiles, notExistsUrls, existsUrls]...
                 =Many( urls, toFolder, dfltFldr, ...
                 waitWhenDone, allowCancel, where)
            if nargin<6
                where='center';
                if nargin<5
                    allowCancel=true;
                    if nargin<4
                        waitWhenDone=false;
                        if nargin<2 || isempty(toFolder)
                            toFolder=WebDownload.DefaultFolder;
                        end
                        if nargin<3
                            dfltFldr=File.TrimHome(toFolder);
                            if strcmp(toFolder, dfltFldr)
                                dfltFldr=WebDownload.DOCUMENTS_FOLDER_FOR_EXAMPLES;
                            end
                            if startsWith(dfltFldr, ['Documents' filesep])
                                dfltFldr=dfltFldr(11:end);
                            end
                        end
                    end
                end
            end
            if ispc % sigh
                dfltFldr=strrep(dfltFldr, filesep, '/');
            end
            File.mkDir(toFolder);
            if ~endsWith(dfltFldr,'/')
                dfltFldr=[dfltFldr '/'];
            end
            ok=false;
            N=length(urls);
            localFiles=cell(1, N);
            for i=1:N
                url=urls{i};
                [p, f, ext]=fileparts(url);
                if isempty(p)
                    url=[WebDownload.HOST dfltFldr f ext];
                    urls{i}=url;
                end
                if startsWith(url, WebDownload.HOST)
                    path=url(length(WebDownload.HOST):end);
                    key=fileparts(path);
                    urls{i}=WebDownload.ResolveUrl([f ext], key);
                end
                localFiles{i}=fullfile(toFolder, ...
                    [String.URLDecode(f) String.URLDecode(ext)]);
            end
            if isempty(localFiles) % nothing to do
                ok=true;
                cancelled=false;
                dwl=[];
                dlg=[];
                notExistsFiles={}; existsFiles={};...
                    notExistsUrls={};existsUrls={};
                return;
            end
            Gui.Modem;
            [cancelled, ~, dwl, dlg]=...
                WebDownload.Get(urls, localFiles,  ...
                waitWhenDone, allowCancel, where);
            if nargout>4
                files = cell(1,N);
                exists = false(1,N);

                for i=1:N
                    [~,f,ext]=fileparts(localFiles{i});
                    files{i} = [f ext];
                    exists(i) = logical(exist(localFiles{i}, 'file'));
                end

                existsFiles = files(exists);
                notExistsFiles = files(~exists);
                existsUrls = urls(exists);
                notExistsUrls = urls(~exists);
            end
            if ~cancelled 
                if ~exist('notExistsUrls', 'var') || isempty(notExistsUrls)
                    ok=true;
                end
            end
        end
        
        %where only used if Gui.m is present
        function [cancelledByUser, bad, dwl, dlg]=...
                Get(urls, localFiles,  waitWhenDone, allowCancel, where)
            if nargin<5
                where='center';
                if nargin<4
                    allowCancel=true;
                    if nargin<3
                        waitWhenDone=true;
                    end
                end
            end
            try
                dwl=javaObjectEDT('edu.stanford.facs.swing.WebDownload');
            catch
                jar=fullfile(fileparts(mfilename('fullpath')), 'webDownload.jar');
                javaaddpath(jar);
                try
                    dwl=javaObjectEDT('edu.stanford.facs.swing.WebDownload');
                catch ex
                    ex.getReport
                    error('Trouble using WebDownload');
                end
            end
            dlg=javaObjectEDT(dwl.dlg);
            try
                Gui.LocateJava(dlg, ...
                    Gui.JWindow(get(0, 'CurrentFigure')), where);
            catch
                Gui.LocateJavaOnScreen(dlg, where);
            end
            dwl.waitWhenDone=waitWhenDone;
            dwl.allowCancel=allowCancel;
            N=length(urls);
            for i=1:N
                if startsWith(urls{i}, 'http://cgworkspace.cytogenie.org')
                    urls{i}=strrep(urls{i}, ...
                        'http://cgworkspace.cytogenie.org', ...
                        'https://storage.googleapis.com/cytogenie.org');
                end
            end
            cancelledByUser=~dwl.go(urls,localFiles, ~isempty(where));
            bad=dwl.bad;
            GoogleDrive.HandleTooBig(dwl, true);
        end
        
        function txt=ReadText(url)
            localFile=[tempname '.txt'];
            [cancelledByUser, bad]=...
                WebDownload.Get({url}, {localFile},  false, false, []);
            if ~cancelledByUser && ~bad
                txt=File.ReadTextFile(localFile);
            else
                txt=[];
            end
        end
        %where only used if Gui.m is present
        function [ok, cancelledByUser]=GetAndUnzip(url, zipFile, waitWhenDone, ...
                allowCancel, where)
            if nargin<5
                where='center';
                if nargin<4
                    allowCancel=true;
                    if nargin<3
                        waitWhenDone=true;
                    end
                end
            end
            ok=false;
            [cancelledByUser, bad, dwl, dlg]=WebDownload.Get({url}, ...
                {zipFile}, waitWhenDone, allowCancel, where);
            if ~cancelledByUser && bad==0
                dlg.setModal(false);
                dwl.progressBar.setValue(0);
                dlg.setVisible(true);
                dwl.south.setText('Unzipping file now');
                try
                    zipFldr=fileparts(zipFile);
                    if isempty(zipFldr)
                        unzip(zipFile);
                    else
                        unzip(zipFile, zipFldr);
                    end
                    dwl.progressBar.setValue(dwl.progressBar.getMaximum);
                    MatBasics.RunLater(@(h,e)quiet,2);
                    ok=true;
                catch ex
                    ex.getReport
                end
            end
            delete(zipFile);
                    
            function quiet
                dlg.setVisible(false);
            end
        end
        
        function [x,y]=LocateJavaOnScreen(javaComponent, where)
           warning('DEPRECATED:  Call same function in Gui.m');
           if nargin<2
               where='center';
           end
           [x,y]=Gui.LocateJavaOnScreen(javaComponent, where);
        end
        
        function [x,y]=LocateWidthHeight(isTop0, width, height, ...
                refX, refY, refWidth, refHeight, where)
           warning('DEPRECATED:  Call same function in Gui.m');
           [x,y]=Gui.LocateWidthHeight(isTop0, width, height, ...
                refX, refY, refWidth, refHeight, where);
        end
        
        function file=GetZipIfMissing(file, url)
            if ischar(file) && ~exist(file, 'file')
                [fldr,fn]=fileparts(file);
                zipFileName=[fn '.zip'];
                if nargin<2
                    url=fullfile(WebDownload.HOST, WebDownload.PATH);
                end
                url=WebDownload.FullUrl(url, zipFileName);
                zipFile=fullfile(fldr, zipFileName);
                [ok, cancelledByUser]=WebDownload.GetAndUnzip(url, ...
                    zipFile, false, true, 'center');
                if cancelledByUser
                     msg(Html.WrapHr(['Cancelling leaves the required '...
                         'file missing...<br>' BasicMap.Global.smallStart...
                         '<b>' file '</b>' BasicMap.Global.smallEnd]));
                     file=[];
                elseif ~ok
                    msg(Html.WrapHr(['Required file is missing...<br>' ...
                        BasicMap.Global.smallStart...
                        '<b>' file '</b>' BasicMap.Global.smallEnd]));
                    file=[];
                end
            end
        end
        
        function url=FullUrl(beforePath, afterPath)
            if endsWith(beforePath, '/')
                beforePath=beforePath(1:end-1);
            end
            if startsWith(afterPath, '/')
                afterPath=afterPath(2:end);
            end
            url=[beforePath '/' afterPath];
        end
        
        function [host, path, port, query, protocol]=UrlParts(url)
            u=matlab.net.URI(url);
            path=char(u.EncodedPath);
            if ispc 
                if ~isempty(path) && path(1)=='/' && ~isempty(find(path==':',1)) %drive letter
                    path=path(2:end);
                end
            end                    
            host=char(u.Host);
            query=char(u.EncodedQuery);
            idx=String.IndexOf(url, ':');
            if idx>0
                protocol=url(1:idx-1);
            else
                protocol='';
            end
            port=u.Port;
        end

        function [fullFile, exists, downloaded, isDifferent]...
                =GetFile(fileName, localFolder, remoteFolder, ...
                force, complainEvenIfExists)
            if nargin<5
                complainEvenIfExists=true;
                if nargin<4
                    force=false;
                end
            end
            fullFile=fullfile(localFolder, fileName);
            File.mkDir(localFolder);
            downloaded=false;
            isDifferent=false;
            exists=exist(fullFile, 'file');
            if force || ~exists
                if nargout>3
                    bak=tempname;
                    copyfile(fullFile, bak);
                end
                url=WebDownload.ResolveUrl(fileName, remoteFolder);
                [cancelled, bad]=WebDownload.Get({url}, {fullFile}, false, false, 'south');
                downloaded= ~cancelled && ~bad;
                if ~downloaded
                    if ~exists || complainEvenIfExists
                        msg(['<html>Could not download the file "' ...
                            fileName '"<br>'...
                            Html.WrapColor('From:  ', 'blue')...
                            Html.WrapSmallTags(url) '<br>' ...
                            Html.WrapColor('To:  ', 'blue') ...
                            Html.WrapSmallTags(fileparts(fullFile)) '<hr></html>']);
                    end
                else
                    exists=true;
                    if nargout>3
                        isDifferent=File.Diff(bak, fullFile);
                    end
                end
            end
        end
        
        function fldr=LocalExamplesFolder
            fldr=fullfile(File.Home, 'Documents', ...
                WebDownload.DOCUMENTS_FOLDER_FOR_EXAMPLES);
            File.mkDir(fldr);
        end
        
        function file=GetExampleIfMissing(file, localExamplesFolder)
            file=File.ExpandHomeSymbol(file);
            if ~exist(file, 'file')
                p=fileparts(file);
                if isempty(p)
                    if nargin<2
                        localExamplesFolder=WebDownload.LocalExamplesFolder;
                    end
                    file=fullfile(localExamplesFolder, file);
                    if ~exist(file, 'file')
                        csvFiles=WebDownload.RelocateExamples(...
                            {file}, true, {}, localExamplesFolder);
                        file=csvFiles{1};
                    end
                end
            end
        end
        
        function [result, existence, missingFiles]=...
                RelocateExamples(files, tryDownload, ignore, localExamplesFolder)
            if nargin<4
                localExamplesFolder=WebDownload.LocalExamplesFolder;
                if nargin<3
                    ignore={};
                    if nargin<2
                        tryDownload=true;
                    end
                end
            end
            testExistence=nargout>1;
            existence=[];
            missingFiles={};
            argType='cell';
            if isstruct(files)
                argType='struct';
                args=files;
                files={};
                if ischar(args.csv_file_or_data)
                    files{end+1}=args.csv_file_or_data;
                end
                if ~isempty(args.label_file)
                    files{end+1}=args.label_file;
                end
                if ~isempty(args.template_file)
                    files{end+1}=args.template_file;
                end
                if ~isempty(args.color_file)
                    files{end+1}=args.color_file;
                end
            elseif ischar(files)
                argType='char';
                files={files};
            end
            missingExamples=java.util.HashMap;
            N=length(files);
            fileUrls=cell(1,N);
            localFiles=cell(1,N);
            nFileUrls = 0;
            for i=1:N
                file=files{i};
                if ~exist(file, 'file')
                    if ~ispc
                        if file(1)=='~'
                            file=[File.Home file(2:end)];
                        end
                    end
                end
                if ~exist(file, 'file')
                    [fldr, fn, ext]=fileparts(file);
                    if isempty(fldr) || isequal(localExamplesFolder, fldr)
                        fldr=localExamplesFolder;
                        File.mkDir(fldr);
                        file2=fullfile(fldr, [fn ext]);
                        missingExamples.put(java.lang.String(file), java.lang.String(file2));
                        
                        if ~exist(file2, 'file') && tryDownload
                            nFileUrls = nFileUrls + 1;
                            fileUrls{nFileUrls}=...
                                WebDownload.ResolveUrl([fn ext]);
                            localFiles{nFileUrls}=file2;
                        end
                    end
                end
            end
            fileUrls = fileUrls(1:nFileUrls);
            localFiles = localFiles(1:nFileUrls);

            downloadFailures=[];
            if ~isempty(fileUrls)
                nMissing=length(fileUrls);
                [cancelledByUser, bad]=WebDownload.Get(...
                    fileUrls, localFiles, false);
                if cancelledByUser
                elseif bad==0
                    msg(Html.WrapHr([ String.Pluralize2('file', ...
                        nMissing) ' downloaded to<br><b>' ...
                        BasicMap.Global.smallStart ...
                        localExamplesFolder...
                        BasicMap.Global.smallEnd '</b>']), 8, ...
                        'south east+');
                else
                    
                    if bad==nMissing
                        downloadFailures=[String.Pluralize2('file', ...
                            nMissing) ' NOT downloaded'];
                    else
                        downloadFailures=[ num2str(bad) ' of '...
                            String.Pluralize2('file', ...
                            nMissing) ' NOT downloaded'];
                    end
                end
            end
            if isequal(argType, 'struct')
                if ischar(args.csv_file_or_data)
                    args.csv_file_or_data=grab(args.csv_file_or_data);
                end
                args.label_file=grab(args.label_file);
                args.template_file=grab(args.template_file);
                args.color_file=grab(args.color_file);
            else
                for i=1:N
                    files{i}=grab(files{i});
                end
            end
            if isequal(argType, 'cell')
                result=files;
            elseif isequal(argType, 'struct')
                result=args;
            else %argType == char
                result=files{1};
            end
            if any(~existence) && ~isempty(missingFiles)
                if ~isempty(ignore)
                    [~, canIgnore]=StringArray.EndsWith(missingFiles, ignore);
                else
                    canIgnore=false;
                end
                if ~canIgnore
                    app=BasicMap.Global;
                    html=Html.Wrap([app.h2Start 'Missing files' ...
                        app.h2End downloadFailures app.smallStart ...
                        Html.ToList(missingFiles, 'ol') app.smallEnd '<hr>']);
                    msgWarning(html, 11, 'south', 'Missing files...');
                else
                    existence(:)=2;
                end
            end
            
            function out=grab(in)
                out=in;
                if ~isempty(in)
                    key=java.lang.String(in);
                    if missingExamples.containsKey(key)
                        if tryDownload
                            out=char(missingExamples.get(key));
                        else
                            out='';
                        end
                    else
                        if in(1)=='~'
                            out=[File.Home in(2:end)];
                        end
                    end
                end
                if testExistence
                    existence(end+1)=exist(out, 'file');
                    if ~existence(end)
                        missingFiles{end+1}=out;
                    end
                end
            end
        end
        
        % convenience method for default local folder
        function fldr=DefaultFolder()
            fldr=fullfile(File.Documents, ...
                WebDownload.DOCUMENTS_FOLDER_FOR_EXAMPLES);
        end

        function [nMissing, missingFiles, columnsAreEmbedded]=GetMlpIfMissing( ...
                modelFile, forPython)
            nMissing=0;
            missingFiles={};
            columnsAreEmbedded=false;
            if isempty(fileparts(modelFile))
                modelFile=fullfile(Mlp.DefaultLocalFolder, modelFile);
            end
            files={[modelFile Mlp.EXT_COLUMNS]};
            if forPython
                files{end+1}=[modelFile Mlp.EXT_TENSORFLOW];
                files{end+1}=[modelFile Mlp.EXT_DICT];
                files{end+1}=[modelFile Mlp.EXT_SCALE];
            else
                files{end+1}=[modelFile Mlp.EXT_FITCNET];
                columnsAreEmbedded=isEmbedded;
                if columnsAreEmbedded
                    return;
                end
            end
            app=BasicMap.Global;
            if app.containsKey(Mlp.PROP_SERVER_FOLDER)
                key=app.get(Mlp.PROP_SERVER_FOLDER);
            else
                key='mlp';
            end
            host=app.get(Mlp.PROP_SERVER_HOST);
            [~, existence, missingFiles]=WebDownload.GetIfMissing( ...
                files, true, {}, key, host);
            if any(~existence)
                key=[];%use folder path from home as key
                [~, existence, missingFiles]=WebDownload.GetIfMissing( ...
                    files, true, {}, key, host);
            end
            nMissing=sum(~existence);
            if nMissing>0
                columnsAreEmbedded=isEmbedded;
                if ~columnsAreEmbedded
                    warning("Mlp model missing %s", ...
                        String.ToString(missingFiles));
                end
            end

            function yes=isEmbedded
                yes=false;
                if ~exist(files{1}, 'file') && exist(files{2})
                    load(files{end}, 'column_names')
                    if exist('column_names', 'var')
                        nMissing=0;
                        missingFiles={};
                        yes=true;
                        return;
                    end
                end
            end
        end

        function [result, existence, missingFiles]=...
                GetIfMissing(files, tryDownload, ...
                ignore, key, host)
            if nargin<5
                host=[];
                if nargin<4
                    key=[];
                    if nargin<3
                        ignore={};
                        if nargin<2
                            tryDownload=true;
                        end
                    end
                end
            end
            testExistence=nargout>1;
            existence=[];
            missingFiles={};
            missingExamples=java.util.HashMap;
            N=length(files);
            fileUrls=cell(1,N);
            localFiles=cell(1,N);
            nFileUrls = 0;
            for i=1:N
                file=File.ExpandHomeSymbol(files{i});
                if isempty(fileparts(file))
                    file=fullfile(WebDownload.LocalExamplesFolder, file);
                    files{i}=file;
                end
                if ~exist(file, 'file')
                    nFileUrls = nFileUrls + 1;
                    [fldr, fn, ext]=fileparts(file);
                    File.mkDir(fldr);
                    missingExamples.put(java.lang.String(file), java.lang.String(file));
                    if isempty(key)
                        if startsWith(fldr, File.Home)
                            key_=fldr(length(File.Home)+2:end);
                        else
                            idx=String.IndexOf(fldr, '/');
                            if idx<1
                                idx=String.IndexOf(fldr, '\');
                            end
                            key_=fldr(idx+1:end);
                        end
                        fileUrls{nFileUrls}=WebDownload.ResolveUrl( ...
                            [fn ext], key_, host);
                    else
                        fileUrls{nFileUrls}=WebDownload.ResolveUrl( ...
                            [fn ext], key, host);
                    end
                    localFiles{nFileUrls}=file;
                end
            end
            fileUrls = fileUrls(1:nFileUrls);
            localFiles = localFiles(1:nFileUrls);
            
            downloadFailures=[];
            if ~isempty(fileUrls)
                nMissing=length(fileUrls);
                [cancelledByUser, bad]=WebDownload.Get(...
                    fileUrls, localFiles, false);
                if cancelledByUser
                elseif bad==0
                    msg(Html.WrapHr([ String.Pluralize2('file', ...
                        nMissing) ' downloaded ']), 8, ...
                        'south east+');
                else
                    if bad==nMissing
                        downloadFailures=[String.Pluralize2('file', ...
                            nMissing) ' NOT downloaded'];
                    else
                        downloadFailures=[ num2str(bad) ' of '...
                            String.Pluralize2('file', ...
                            nMissing) ' NOT downloaded'];
                    end
                end
            end
            for i=1:N
                files{i}=grab(files{i});
            end
            result=files;
            if any(~existence) && ~isempty(missingFiles)
                if ~isempty(ignore)
                    [~, canIgnore]=StringArray.EndsWith(missingFiles, ignore);
                else
                    canIgnore=false;
                end
                if ~canIgnore
                    app=BasicMap.Global;
                    html=Html.Wrap([app.h2Start 'Missing files' ...
                        app.h2End downloadFailures app.smallStart ...
                        Html.ToList(missingFiles, 'ol') app.smallEnd '<hr>']);
                    msgWarning(html, 11, 'south', 'Missing files...');
                else
                    existence(:)=2;
                end
            end
            
            function out=grab(in)
                out=in;
                if ~isempty(in)
                    fileKey=java.lang.String(in);
                    if missingExamples.containsKey(fileKey)
                        if tryDownload
                            out=char(missingExamples.get(fileKey));
                        else
                            out='';
                        end
                    else
                        if in(1)=='~'
                            out=[File.Home in(2:end)];
                        end
                    end
                end
                if testExistence
                    existence(end+1)=exist(out, 'file');
                    if ~existence(end)
                        missingFiles{end+1}=out;
                    end
                end
            end
        end
        
        function file=LocateUri(uri, downloadsFolder, fromWspFile, ask)
            if nargin<4
                ask=true;
                if nargin<3
                    fromWspFile=false;
                    if nargin<2
                        downloadsFolder=[];
                    end
                end
            end
            u=lower(uri);
            download=false;
            decodeFileUrl=false;
            if startsWith(u, 'file://')
                file=uri(8:end);
                decodeFileUrl=true;
             elseif startsWith(u, 'file:/')
                file=uri(7:end);
                decodeFileUrl=true;
            elseif startsWith(u, 'file:')
                file=uri(6:end);
                decodeFileUrl=true;
            elseif startsWith(u, 'https://') 
                download=true;
                file=uri(9:end);
            elseif startsWith(u, 'http://')
                download=true;
                file=uri(8:end);
            else
                file=uri;
            end
            if download
                [host, p, f, e]=WebDownload.CloudUriToFile( ...
                    uri, fromWspFile);
                if isempty(downloadsFolder)
                    p=strrep(p, '.app/', '_app/');
                    p=File.Downloads('suh_pipelines', host, p);
                    if exist(fullfile(p, [f e]), 'file')
                        file=fullfile(p, [f e]);
                        return;
                    end
                    if exist(p, 'dir')
                        file=fullfile(p, [f e]);
                    else
                        File.mkDir(p);
                        if ~ask
                            folder=p;
                        else
                            folder=uiPutFile(p, f, BasicMap.Global, ...
                                'SuhPipelines.DownLoads', ...
                                ['Pick "' host '" downloads folder'],...
                                false);
                        end
                        if ~isempty(folder)
                            file=fullfile(folder, [f e]);
                            if exist(file, 'file')
                                [yes, cancelled]=askYesOrNo(...
                                    'Re-download the file?', ...
                                    'File exists ...');
                                if cancelled
                                    file=[];
                                    return;
                                end
                                if ~yes
                                    return;
                                end
                            end
                        else
                            file=[];
                            File.rmDir(p);
                            return
                        end
                    end
                else
                    file=fullfile(downloadsFolder, [f e]);
                    if exist(file, 'file') 
                        return;
                    end
                end
                if fromWspFile &&  endsWith(lower(file), '.wsp')
                    uri2=[uri(1:end-4) '.guiV2.properties'];
                    file2=[file(1:end-4) '.guiV2.properties'];
                    [p,f,e]=fileparts(file2);
                    cancelled=WebDownload.Get({uri, uri2}, {file, file2}, false);
                    if cancelled || ~exist(file, 'file')
                        file=[];
                    end
                    if exist(file2, 'file')
                        try
                            % in case it is not linked
                            fldr=FlowJoWsp.SUH_FOLDER;
                        catch
                            fldr='FlowJoBridgeHerzenbergLab';
                        end
                        [~,f]=fileparts(f);
                        fldr=fullfile(p,  f, fldr);
                        if ~exist(fldr,'dir')
                            File.mkDir(fldr);
                        end
                        file3=fullfile(fldr, 'gui.properties');
                        movefile(file2, file3, 'f');
                    end
                else
                    [cancelled, bad]=WebDownload.Get({uri}, {file}, false);
                    if cancelled || bad
                        file=[];
                    end
                end                
            elseif decodeFileUrl
                file=WebDownload.FileUriToFile(uri, fromWspFile);
            end
        end
    end
    
    properties(Constant)
        FLOWJO_FILE_UNENCODED='&+;~=!@()';
        %FLOWJO_FILE_ENCODED=urlencode(WebDownload.FLOWJO_FILE_UNENCODED);
        CACHE_DIR='/Downloads/suh_pipelines/';
        GOOGLE_CYTOGENIE=[WebDownload.CACHE_DIR 'storage.googleapis.com/'];
    end

    methods(Static)
        function [host, path, fileName, extension]...
            =CloudUriToFile(uri, fromWspFile)
            if fromWspFile
                [~, host, path]=WebDownload.FlowJoFileUriFromCloud(uri);
                path=WebDownload.FileUriToFile(path, true);
            else
                [host, path]=WebDownload.GetHostPath(uri);
                path=WebDownload.UrlDecode(path);
            end
            [path, fileName, extension]=fileparts(path);
        end
        
        function file=FileUriToFile(uri, fromFlowJoWsp)
            if startsWith(uri, 'file:') && contains(uri, ...
                    [WebDownload.GOOGLE_CYTOGENIE 'cytogenie.app'])
                %support unfortunate old BSTR style hack or PSP style hack
                warning('Supporting early hack in cytogenie.app folder caching with %s ', uri);
                uri=strrep(uri, 'cytogenie.app', 'cytogenie_app');%sigh
            end
            if contains(uri, WebDownload.CACHE_DIR)
                idx=String.IndexOf(uri, WebDownload.CACHE_DIR);
                if ~contains(uri, strrep(File.Home, '\', '/'))
                    warning('RELOCATING %s ', uri);
                end
                uri=fullfile(File.Home, uri(idx:end));
            elseif startsWith(lower(uri), 'file://')
                uri=uri(8:end);
            elseif startsWith(lower(uri), 'file:/') ...
                    && uri(8)==':' %stupid Microsoft drive letter!
                if ispc
                    uri=uri(7:end);
                else 
                    uri=uri(9:end);
                end
            elseif startsWith(lower(uri), 'file:')
                uri=uri(6:end);
            end
            if fromFlowJoWsp
                BAD=WebDownload.FLOWJO_FILE_UNENCODED;
                N=length(BAD);
                for i=1:N
                    bad=BAD(i);
                    good=WebDownload.UrlEncode(bad);
                    uri=strrep(uri, bad, good);
                end
            end
            file=WebDownload.UrlDecode(uri);
             if ispc
                file=strrep(file, '/', '\');
            end
        end

        function [host, path]=GetHostPath(uri)
            u=uri(String.IndexOf(uri, '://')+3:end);
            idx=String.IndexOf(u, '/');
            host=u(1:idx-1);
            path=u(idx:end);
        end
        
        function [fileUri, host, path]=FlowJoFileUriFromCloud(uriCloud, rootFolder)
            orig=uriCloud;
            if nargin<2
                rootFolder=File.Downloads('suh_pipelines');
            end
            BAD=WebDownload.FLOWJO_FILE_UNENCODED;
            N=length(BAD);
            for i=1:N
                bad=BAD(i);
                good=WebDownload.UrlEncode(bad);
                uriCloud=strrep(uriCloud, good, bad);
            end
            [host,path]=WebDownload.GetHostPath(uriCloud);
            if ispc
                rootFolder=strrep(rootFolder, '\', '/');
            end
            fileUri=['file:' fullfile(rootFolder, host, path)];
            if ispc
                fileUri=strrep(fileUri, '\', '/');
            end
            if contains(fileUri, ...
                    ['Downloads/suh_pipelines/'...
                    'storage.googleapis.com/cytogenie.app'])
                %support unfortunate old BSTR style hack or PSP style hack
                fileUri=strrep(fileUri, 'cytogenie.app', 'cytogenie_app');%sigh
            end
            
        end

        function cloudUri=FlowJoFileUriToCloud(uriFile)
            if ispc && startsWith(uriFile, 'file:/')
                rootFolder=['file:/' File.Downloads('suh_pipelines', ...
                    'storage.googleapis.com')];
            else
                rootFolder=['file:' File.Downloads('suh_pipelines', ...
                    'storage.googleapis.com')];
            end
            if ispc
                rootFolder=strrep(rootFolder, '\', '/');
            end
            uriFile=uriFile(length(rootFolder)+2:end);
            BAD=WebDownload.FLOWJO_FILE_UNENCODED;
            N=length(BAD);
            for i=1:N
                bad=BAD(i);
                good=WebDownload.UrlEncode(bad);
                uriFile=strrep(uriFile, bad, good);
            end
            cloudUri=['https://storage.googleapis.com/' uriFile];
        end

        function urlOut = UrlEncode(urlIn)
            urlOut = char(java.net.URLEncoder.encode(urlIn,'UTF-8'));
        end

        function urlOut = UrlDecode(urlIn)
            urlOut = char(java.net.URLDecoder.decode(urlIn,'UTF-8'));
        end
    end
    
end