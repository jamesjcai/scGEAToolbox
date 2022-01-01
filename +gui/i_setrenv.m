function i_setrenv(~,~)

%see also: I_SETPYENV
Rpath=pkg.FindRpath;
if isempty(Rpath)
    warndlg('R is not installed.','R Environment');
else
    if iscell(Rpath)
        s=Rpath{1};
        for k=2:length(Rpath)
            s=sprintf('%s\n%s',s,Rpath{k});
        end
        s=sprintf('%s %s',s,' (default)');
    else
        s=Rpath;
    end
    helpdlg(sprintf('R executable found at:\n%s',...
            s),'R Environment');
end
end