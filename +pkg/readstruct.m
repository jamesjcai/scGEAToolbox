function S = readstruct(filename,varargin)
%READSTRUCT Create a struct by reading from a file.
%
%   Use the READSTRUCT function to create a struct by reading
%   structured data from a file. READSTRUCT automatically determines
%   the file format from its extension.
%
%   S = READSTRUCT(FILENAME) creates a struct by reading from a file, where
%   FILENAME can be one of these:
%
%       - For local files, FILENAME can be an absolute path that contains
%         a filename and file extension. FILENAME can also be a relative path
%         to the current folder, or to a folder on the MATLAB path.
%         For example, to import a file on the MATLAB path:
%
%            S = readstruct("music.xml");
%
%       - For files from an Internet URL or stored at a remote location,
%         FILENAME must be a full path using a Uniform Resource Locator
%         (URL). For example, to import a remote file from Amazon S3,
%         specify the full URL for the file:
%
%            S = readstruct("s3://bucketname/path_to_file/my_setup.xml");
%
%         For more information on accessing remote data, see "Work with
%         Remote Data" in the documentation.
%
%   S = READSTRUCT(FILENAME,"FileType",FILETYPE) specifies the file type, where
%   the default FILETYPE is "auto". This only needs to be specified if FILENAME
%   does not have a recognized extension (e.g. an XML file with a .gpx extension).
%
%   Name-Value Pairs:
%   -------------------------------------------------------------------------------
%
%   "FileType"             - Specifies the file type. It can be specified as:
%                            - "auto": infers the FileType from the file extension.
%                            - "xml":  treat as an XML file, regardless of file
%                                      extension.
%
%   "StructNodeName"       - Name of XML Element node underneath which READSTRUCT
%                            should start reading a struct.
%
%   "StructSelector"       - XPath expression which precisely specifies the XML
%                            Element node underneath which READSTRUCT should start
%                            reading in a struct.
%
%   "ImportAttributes"     - Import XML node attributes as fields of the output
%                            struct. Defaults to true.
%
%   "AttributeSuffix"      - Suffix to append to all output struct field names
%                            corresponding to attributes in the XML file. Defaults
%                            to "Attribute".
%
%   "RegisteredNamespaces" - The namespace prefixes that are mapped to
%                            namespace URLs for use in selector expressions.
%
%   "DateLocale"           - The locale used to interpret month and day
%                            names in datetime text. Must be a character
%                            vector or a scalar string in the form xx_YY.
%                            See the documentation for DATETIME for more
%                            information.
%
%   "WebOptions"           - HTTP(s) request options, specified as a 
%                            weboptions object. 
%
%   See also WRITESTRUCT, STRUCT

%   Copyright 2019-2022 The MathWorks, Inc.

    try
        func = matlab.io.internal.functions.FunctionStore.getFunctionByName("readstruct");
        S = func.validateAndExecute(filename,varargin{:});
    catch ME
        throw(ME);
    end
end
