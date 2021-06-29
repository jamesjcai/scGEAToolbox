function i_setrenv(~,~)
    helpdlg(sprintf('R executable at:\n%s',...
        FindRpath),'R Environment');
end