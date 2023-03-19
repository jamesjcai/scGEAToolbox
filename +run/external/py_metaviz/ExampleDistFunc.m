%This is nonsense distance provided here to illustrate callback API
function D2=ExampleDistFunc(Z1,ZJ)
%     [R,C]=size(ZJ);
%     for r=1:R
%         if all(ZJ(r,:)==Z1)
%             D2=ZJ(r,:)';
%         end
%     end
    z1=Z1(1:50);
    zj=ZJ(:,1:50);
    [~, index]=ismember(z1,zj,'rows');
    D2=ZJ(:,index);
end
