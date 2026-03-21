function [a] = i_guessspecies(g)

a = 'human';
if sum(upper(g) == g)./numel(g) < 0.9
    a = 'mouse';
end
