function i_setrenv(~,~)
    helpdlg(sprintf('R executable at:\n%s',...
        pkg.FindRpath),'R Environment');
end