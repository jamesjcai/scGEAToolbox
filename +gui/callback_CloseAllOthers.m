function callback_CloseAllOthers(~,~)
a=gcf;
all_figs = findobj(0, 'type', 'figure');
[y,idx]=ismember(a,all_figs);
if y && length(all_figs)>1
   answer=questdlg('Close all other figures?');
   if ~strcmp(answer,'Yes'), return; end
   for k=1:length(all_figs)
       if k~=idx
           try
            close(all_figs(k));
           catch
           end
       end
   end
end
end