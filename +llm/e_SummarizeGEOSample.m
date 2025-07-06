prompt = "Summarize the sample found at the provided link https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7855468 from NCBI GEO. Return only the summary text.";

% response = llm.geminiGenerateContent(prompt);
% if response.StatusCode == "OK"
%     response.Body.Data.candidates.content.parts.text
% else
%     response.Body.Data.error
% end

    preftagname = 'llmodelprovider';
    s = getpref('scgeatoolbox', preftagname);
    providermodel = strsplit(s,':');

    provider = providermodel{1};
    fprintf('Using LLM provider: %s\n', provider);
    model = providermodel{2};
    fprintf('Using LLM model: %s\n', model);


[done, response] = llm.callTAMUAIChat([], prompt, model)
