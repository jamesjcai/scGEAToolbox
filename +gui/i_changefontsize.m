function i_changefontsize(~, ~)
ax = get(gca, 'FontSize') - 1;
if ax <= 5, ax = 15; end
set(gca, 'FontSize', ax);
%      ax=get(gca,'LabelFontSizeMultiplier')/1.1;
%      if ax<=0.5, ax=3; end
%      set(gca,'LabelFontSizeMultiplier')
end
