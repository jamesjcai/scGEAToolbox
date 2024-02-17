function [str] = num2bankScalar(num)
    % https://www.mathworks.com/matlabcentral/answers/96131-is-there-a-format-in-matlab-to-display-numbers-such-that-commas-are-automatically-inserted-into-the
    num = floor(num*100) / 100;
    str = num2str(num);
    k = find(str == '.', 1);
    if isempty(k)
        % str=[str,'.00'];
    end
    % FIN = min(length(str),find(str == '.')-1);
    FIN = length(str);
    for i = FIN - 2:-3:2
        str(i + 1:end + 1) = str(i:end);
        str(i) = ',';
    end
end
