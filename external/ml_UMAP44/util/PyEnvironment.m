classdef PyEnvironment < handle
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause

    methods(Static)
        function ok=Setup
            if ~verLessThan('matLab', '9.8')
                app=BasicMap.Global;
                if ~app.pyEnvRun
                    try
                        jd=msg('Initializing Python TensorFlow', 0);
                        app.pyEnvRun=true;
                        directory = fullfile(fileparts( ...
                            mfilename('fullpath')), 'mlp');
                        P = py.sys.path;
                        if count(P,directory) == 0
                            insert(P,int32(0),directory);
                        end
                        py.importlib.import_module('mlp');
                        ok=true;
                    catch ex
                        if ~isdeployed
                            ex.getReport
                        end
                    end
                    if exist('jd', 'var')
                        jd.dispose;
                    end
                end
            end
        end
    end
end