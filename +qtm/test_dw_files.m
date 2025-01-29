% chat = ollamaChat("deepseek-r1", TimeOut = 1200);
chat = ollamaChat("qwen2.5-coder", TimeOut = 1200);

url = 'https://github.com/jamesjcai/scGEAToolbox/tree/main/%2Bgui';
html = webread(url);

% Extract file links (simplified approach, may need regex adjustments)
file_links = regexp(html, 'href="(/jamesjcai/scGEAToolbox/blob/main/%2Bgui/.*?\.m)"', 'tokens');
file_links = [file_links{:}];
file_links = unique(file_links);

base_url = 'https://raw.githubusercontent.com/jamesjcai/scGEAToolbox/main/%2Bgui/';

for i = 1:2 %length(file_links)
    filename = extractAfter(file_links{i}, '%2Bgui/')
    file_url = [base_url filename];
    outfile = filename(1:end-2) + "_coder.txt";
    fid = fopen(outfile, 'w');
    try
        code = webread(file_url);
        prompt = "Review this MATLAB code for professional feedback: " + code;
        feedbk = generate(chat, prompt);        
        fprintf(fid, '%s', feedbk);
    catch
        fprintf('Failed to read %s\n', filename);
    end
    fclose(fid);
end

%{
chat = ollamaChat("deepseek-r1", StreamFun=@printToken);
prompt = "What is Model-Based Design?";
generate(chat, prompt, MaxNumTokens=500);
function printToken(token)
    fprintf("%s",token);
end
%}