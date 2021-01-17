function i_violinplot_groupordered2(d,Gender,Smoker)
import pkg.Violin
import pkg.violinplot

if ~isstring(Gender)
    Gender=string(Gender);
end
if ~isstring(Smoker)
    Smoker=string(Smoker);
end

[G,gender,smoker] = findgroups(Gender,Smoker);
%meanSystolic = splitapply(@mean,Systolic,G);
%meanDiastolic = splitapply(@mean,Diastolic,G);

Gender=strrep(Gender,'_',' ');
[~,cL]=grp2idx(Gender);
[~,i]=sort(grpstats(d,Gender,@median),'descend');

violinplot(d,Gender,'GroupOrder',cL(i),...
    'ShowData',false,... %'ViolinColor',[1 1 1],...
    'EdgeColor',[0 0 0]);
xtickangle(-45);
box on
ylabel('Differentiation Potency');
grid on