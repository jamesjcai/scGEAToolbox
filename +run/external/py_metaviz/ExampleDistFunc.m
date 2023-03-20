function D2=ExampleDistFunc(Z1,ZJ)
    %     [R,C]=size(ZJ);
    %     for r=1:R
    %         if all(ZJ(r,:)==Z1)
    %             D2=ZJ(r,:)';
    %         end
    %     end
    n=min([50 length(Z1)]);
    z1=Z1(1:n);
    zj=ZJ(:,1:n);
    [~, index]=ismember(z1,zj,'rows');
    D2=ZJ(:,index);
end
