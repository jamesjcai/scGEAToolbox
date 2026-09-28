function myWarndlg(parentfig, message, title, modal)

if nargin < 4, modal = true; end
if nargin < 3, title = ''; end
if nargin < 2, message = 'Message'; end
if nargin < 1, parentfig = []; end
if gui.i_isuifig(parentfig)
    tag = 'warning';
else
    tag = 'warn';
end
gui.mydlg(tag, parentfig, message, title, modal);

