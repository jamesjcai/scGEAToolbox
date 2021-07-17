function i_setrenv(~,~)

%see also: I_SETPYENV
if isempty(pkg.FindRpath)
    warndlg('R is not installed.','R Environment');
else
    helpdlg(sprintf('R executable found at:\n%s',...
        pkg.FindRpath),'R Environment');
end
end