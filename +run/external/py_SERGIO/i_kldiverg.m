function [kl]=i_kldiverg(pt,po)
    KL1 = sum(po .* (log(po)-log(pt)));
    KL2 = sum(pt .* (log(pt)-log(po)));
    kl=10*mean([KL1 KL2]);
end