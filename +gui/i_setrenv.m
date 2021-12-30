function i_setrenv(~,~)

%see also: I_SETPYENV
Rpath=pkg.FindRpath;
if isempty(Rpath)
    warndlg('R is not installed.','R Environment');
else
    if iscell(Rpath)
        Rpath=Rpath{end};
    end 
    helpdlg(sprintf('R executable found at:\n%s',...
        Rpath),'R Environment');
end
end