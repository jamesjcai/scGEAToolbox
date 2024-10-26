function out=uiGetFile(clue, folder, ttl, properties, property)
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause

out=[];
if nargin<2 || isempty(folder)
    folder=File.Documents;
end
if nargin>3
    fldr=properties.get(property, folder);
    if ~isempty(fileparts(fldr))
        folder=fldr;
    end
end
if ismac
    jd=Gui.MsgAtTopScreen(ttl,25);
else
    jd=[];
end
if exist(folder, 'dir')
    [file, fldr]=uigetfile(clue, char(...
        edu.stanford.facs.swing.Basics.RemoveXml(ttl)), ...
        [folder '/']);
elseif exist(folder, 'file')
    [file, fldr]=uigetfile(clue, char(...
        edu.stanford.facs.swing.Basics.RemoveXml(ttl)), ...
        folder);
else
    [file, fldr]=uigetfile(clue, char(...
        edu.stanford.facs.swing.Basics.RemoveXml(ttl)), ...
        fileparts(folder));
end
if ~isempty(jd)
    jd.dispose;
end
if ~isnumeric(file) && ~isnumeric(fldr)
    out=fullfile(fldr,file);
end
if isempty(out)
    return;
end
if nargin>3
    properties.set(property, out);
end

end
        