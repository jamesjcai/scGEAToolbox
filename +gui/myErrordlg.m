function [h] = myErrordlg(parentfig, message, title, modal)

if nargin < 4, modal = true; end
if nargin < 3, title = ''; end
if nargin < 2, message = 'Message'; end
if nargin < 1, parentfig = []; end
gui.mydlg('error', parentfig, message, title, modal);

