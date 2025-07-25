function [folder, file]=uiPutFile(proposedFldr, ...
    proposedFile, propertiesObject, propertyNames, theTitle, ...
    proposalIsDefault, filters, reOrganizeFilters)
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause

if nargin<8
    reOrganizeFilters=false;
    if nargin<7
        filters={};
        if nargin<6
            proposalIsDefault=true;
            if nargin<5
                theTitle='Save to which folder & file?';
                if nargin<4
                    propertyNames='uiPutFile';
                    if nargin<3
                        propertiesObject=BasicMap.Global;
                    end
                end

            end
        end
    end
end
if isempty(proposedFldr)
    proposedFldr=File.Documents;
end
File.mkDir(proposedFldr);
propFile='';
if iscell(propertyNames)
    propFolder=propertyNames{1};
    propFile=propertyNames{2};
else
    propFolder=propertyNames;
end
if ~isempty(propertiesObject)
    lastFolder=propertiesObject.get(propFolder, proposedFldr);
    if ~exist(lastFolder, 'dir')
        if exist(fileparts(lastFolder), 'file')
            lastFolder=fileparts(lastFolder);
        else
            lastFolder = proposedFldr;
        end
    end
else
    lastFolder=proposedFldr;
end
startingFldr=lastFolder;
if proposalIsDefault
    dfltFldr=proposedFldr;
else
    dfltFldr=lastFolder;
end
if ~isempty(propFile) && ~isempty(propertiesObject)
    proposedFile=propertiesObject.get(propFile, proposedFile);
end
[~,~,ext]=fileparts(proposedFile);
done=false;
jd=Gui.MsgAtTopScreen(theTitle, 25);
if startsWith(theTitle, '<html>')
    theTitle=char(edu.stanford.facs.swing.Basics.RemoveXml(theTitle));
end
if reOrganizeFilters 
    filters=File.ShuffleFilters(filters, ext);
end
while ~done
    done=true;
    if isempty(filters)
        [file, folder]=uiputfile(['*' ext], ...
            theTitle, fullfile(startingFldr, proposedFile));
    else
        [file, folder]=uiputfile(filters, ...
            theTitle, fullfile(startingFldr, proposedFile));
    end
    if ~isempty(jd)
        jd.dispose;
    end
    if isempty(folder) || isnumeric(folder)
        folder=[];
        file=[];
        if isequal(dfltFldr, startingFldr)
            return;
        end
        if isequal([dfltFldr filesep], startingFldr)
            return;
        end
        if isequal(dfltFldr, [startingFldr filesep])
            return;
        end
        if ~File.WantsDefaultFolder(dfltFldr) 
            return;
        end
        [file, folder]=uiputfile(['*' ext], ...
            'Save to which folder & file?', ...
            fullfile(proposedFldr, proposedFile));
        if isempty(folder)|| isnumeric(folder)
            folder=[];
            file=[];
            return;
        end
    end
end
propertiesObject.set(propFolder, folder);
if ~isempty(propFile)
    propertiesObject.set(propFile, file);
end

    
end