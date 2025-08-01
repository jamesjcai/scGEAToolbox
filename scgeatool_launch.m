function varargout = scgeatool_launch(varargin)
% SCGEATOOL - Launch either original figure GUI or App Designer version
%
% Use:
%   scgeatool            % launches  App Designer by default
%   scgeatool('-legacy') % launches older figure-based GUI
%   scgeatool(..., 'foo', bar, ...) applies inputs to chosen version

    % Look for 'Legacy' flag
    isLegacy = any(strcmpi(varargin, '-legacy'));

    if isLegacy
        % Remove the flag before calling legacy function
        args = varargin(~strcmpi(varargin, '-legacy'));
        f = scgeatool_legacy(args{:});
        if nargout > 0
            varargout{1} = f;
        end
    else
        try
            % Direct constructor call to App Designer app
            app = scgeatoolApp(varargin{:});
            if nargout > 0
                varargout{1} = app;
            end
        catch ME
            % Fallback if direct call fails
            try
                run('scgeatoolApp.mlapp');
            catch
                rethrow(ME);
            end
        end
    end
end
