%SHAREDDATA enables sharing of a 4x4 transform across processes.
%   SHAREDDATA(name, class, size) provides access to data that can be
%   shared across processes, e.g., across multiple MATLAB instances. This
%   is done using memory mapped files, which results in practically zero
%   overhead.
%   
%   Notes:
% 
%   -) Memory mapping requires a binary file backing the accessed memory.
%   The backing file is called '<name>.tmp' and is located in the directory
%   of this class. The backing file gets created by the first caller. The
%   file is deleted automatically once there are no more callers.
%
%   -) The data gets initialized to all zeros on file creation
%
%   -) There can be multiple reader threads (and processes), but there
%   should only ever be a single writer thread.
%
%   -) Reads between threads are not guaranteed as MATLAB does not provide
%   access to memory barriers. Reads may also interfere with writes, as
%   there is no locking.
%
%   -) The transform can be shared with processes in other languages. 
%   However, keep in mind that MATLAB stores matrices in column major 
%   format with potentially different endian-ness.
%
%   Example:
%      % MATLAB instance (writer)
%      data = [1 2 3];
%      output = SharedData('vec3', class(data), size(data));
%      output.data = data;
%
%      % MATLAB instance (reader)
%      input = SharedData('vec3', 'double', [1 3]);
%      data = input.data;
%

%   Author: Florian Enner <florian @ hebirobotics.com>
%   Copyright 2015 HEBI Robotics, LLC.
classdef SharedData < handle
    
    properties (SetAccess = private)
        store % path to the backing file
        bytes
        size
        class
    end
    
    properties (Access = public)
        data
    end
    
    properties (Access = private)
        file
        fileLock
    end
    
    methods
        
        % constructor
        function this = SharedData(name, precision, dim)
            
            this.class = precision;
            this.size = dim;
            emptyData = zeros(dim, precision);
            w = whos('emptyData');
            this.bytes = w.bytes;
            
            % find where this file is located
            scriptPath = which(mfilename());
            [path,~,~] = fileparts(scriptPath);
            
            % find the corresponding temporary file
            % in the directory of this script
            fileName = fullfile(path, [name '.tmp']);
            
            % make sure the file exists (create on first call)
            % (2 is magic number for file, not folder)
            if(exist(fileName, 'file') ~= 2)
                display(['creating backing file ' fileName]);
                
                % create file initialized to zeros
                fileID = fopen(fileName,'w');
                fwrite(fileID, emptyData, this.class);
                fclose(fileID);
            end
            
            % memmap does not keep the file from being deleted, so we 
            % need to open the file manually in order to provide
            % auto-deletion.
            this.store = fileName;
            this.fileLock = fopen(this.store,'r');
            
            % make sure the fileSize is correct
            fileInfo = dir(fileName);
            if(fileInfo.bytes ~= this.bytes)
                error(['Backing file has wrong dimensions. '...
                    'Please delete: ' fileName]);
            end
            
            % memory map the file
            this.file = memmapfile(fileName, ...
                'Format', {this.class this.size 'data';}, ...
                'Writable', true);
            
        end
        
        % destructor
        function delete(this)
            
            % unmap memory
            this.file = [];
            
            % give up lock and try to delete file. Note that we disable
            % warnings, as only the last one has permission to actually
            % delete the file.
            fclose(this.fileLock);
            warning('off', 'MATLAB:DELETE:Permission');
            delete(this.store);
            warning('on', 'MATLAB:DELETE:Permission');
            
        end
        
        function set.data(this, value)
            if( ~isequal(size(value), this.size) || ~isa(value, this.class))
               error('invalid size or class'); 
            end
            this.file.Data.data = value * 1;
        end
        
        function [data] = get.data(this)
            data = this.file.Data.data * 1;
        end
        
    end
end
