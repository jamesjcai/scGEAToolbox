function response = geminiGenerateContent(prompt, nvp)
    arguments
        prompt      (1,1) {mustBeTextScalar,mustBeNonempty}
        nvp.image   (1,1) {mustBeTextScalar} = "";
    end

    if strlength(nvp.image) == 0
        model = "gemini-2.0-flash";
        query = struct("contents",[]);
        query.contents = {struct("parts",[])};
        query.contents{1}.parts{1} = {struct("text",prompt)};
    else
        model = "gemini-pro-vision";
        if contains(nvp.image,["http://","https://"])
            imdata = imread(nvp.image);
            imwrite(imdata,"imdata.png")
            img = "imdata.png";
        else
            img = nvp.image;
        end
        fid = fopen(img);
        im = fread(fid,'*uint8');
        fclose(fid);
        b64 = matlab.net.base64encode(im);
        [~,~,ext] = fileparts(img);
        MIMEType = "image/" + erase(ext,".");
        query = struct("contents",[]);
        query.contents = {struct("parts",[])};
        query.contents{1}.parts = {struct("text",prompt),struct("inline_data",[])};
        query.contents{1}.parts{2}.inline_data = struct("mime_type",MIMEType,"data",[]);
        query.contents{1}.parts{2}.inline_data.data = b64;
        if isfile("imdata.png")
            delete("imdata.png")
        end
    end
    %https://generativelanguage.googleapis.com/v1beta/models/gemini-2.0-flash:generateContent?key=

    endpoint = "https://generativelanguage.googleapis.com/v1beta/";
    method = "generateContent";
    
    import matlab.net.*
    import matlab.net.http.*  
    apikey = getenv("GEMINI_API_KEY");
    headers = HeaderField('Content-Type', 'application/json');
    request = RequestMessage('post', headers, query);
    response = send(request, URI(endpoint + ...
        "models/" + model + ...
        ":" + method + ...
        "?key=" + apikey));
end