% loadenv("gemini_api_key.env")
% prompt = "Tell me 5 jokes";
prompt = "Why is the sky blue?";
response = llm.geminiGenerateContent(prompt);
if response.StatusCode == "OK"
    response.Body.Data.candidates.content.parts.text
else
    response.Body.Data.error
end