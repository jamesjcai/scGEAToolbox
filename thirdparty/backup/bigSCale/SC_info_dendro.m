function ix=SC_info_dendro( outperm,indici )


for k=1:length(indici)
  start(k) = min(find(ismember(outperm,indici{k})));
end

[c ix]=sort(start);

%disp('Ordine dei clusters:');
%ix

end

