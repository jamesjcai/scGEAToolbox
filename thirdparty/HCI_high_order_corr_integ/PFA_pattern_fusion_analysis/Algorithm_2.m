function w_1= Algorithm_2( sample_num, w, E_1, E_2, E_3, lam_1)
%% fix Y, solve W based Algorithm 2
E_or = [E_1;E_2;E_3];
E = (E_or)/sum(E_or);
[E_add,index]=sort(E);

p=inf;

for i=3*sample_num:-1:1
    o=(2*lam_1+sum(E_add(1:i,1)))/i - E_add(i);
    
    if o>=0
        p=i;
        break
    end
end

o=(2*lam_1+sum(E_add(1:p,1)))/p;
w(1:p,1)=(o-E_add(1:p,1))/(2*lam_1);
w((p+1):(3*sample_num))=0;

w_1 = -100*zeros(size(w,1),1);

w_1(index)=w;

end

