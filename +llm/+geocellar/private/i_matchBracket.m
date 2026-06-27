function endPos = i_matchBracket(text, startPos, openCh, closeCh)
depth = 0;
inStr = false;
esc   = false;
endPos = -1;
for k = startPos:numel(text)
    c = text(k);
    if esc
        esc = false;
        continue;
    end
    if inStr
        if c == '\'
            esc = true;
        elseif c == '"'
            inStr = false;
        end
    else
        if c == '"'
            inStr = true;
        elseif c == openCh
            depth = depth + 1;
        elseif c == closeCh
            depth = depth - 1;
            if depth == 0
                endPos = k;
                return;
            end
        end
    end
end
end
