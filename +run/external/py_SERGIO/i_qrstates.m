function [p]=i_qrstates(C,theta)
    C0=C;
    for k=(4+1):length(C0.Gates)
        C0.Gates(k).Angles=theta(k-4);
    end
    
    S0 = simulate(C0);
    [~,p0] = querystates(S0);
    
    S = simulate(C);
    [~,p] = querystates(S);    
    [p0 p]'
end