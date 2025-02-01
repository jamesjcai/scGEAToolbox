% chat = ollamaChat("deepseek-r1", TimeOut = 1200);
chat = ollamaChat("qwen2.5-coder", TimeOut = 1200);

url = 'https://github.com/jamesjcai/scGEAToolbox/tree/main/%2Bgui';
html = webread(url);

% Extract file links (simplified approach, may need regex adjustments)
file_links = regexp(html, 'href="(/jamesjcai/scGEAToolbox/blob/main/%2Bgui/.*?\.m)"', 'tokens');
file_links = [file_links{:}];
file_links = unique(file_links);

base_url = 'https://raw.githubusercontent.com/jamesjcai/scGEAToolbox/main/%2Bgui/';

N = length(file_links);
fw = gui.gui_waitbar_adv;
for i = 51:N
    gui.gui_waitbar_adv(fw, i/N);
    filename = extractAfter(file_links{i}, '%2Bgui/');
    file_url = [base_url filename];
    outfile = filename(1:end-2) + "_coder.txt";
    if ~exist(outfile, 'file')
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
    
end
gui.gui_waitbar_adv(fw);
%{
chat = ollamaChat("deepseek-r1", StreamFun=@printToken);
prompt = "What is Model-Based Design?";
generate(chat, prompt, MaxNumTokens=500);
function printToken(token)
    fprintf("%s",token);
end
%}