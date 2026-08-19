function h = i_toHash(txt)
% --- Helper: SHA-256 Hashing (no Java runtime required) ---
digester = matlab.internal.crypto.BasicDigester("SHA256");
d = digester.computeDigest(uint8(txt));
h = lower(dec2hex(d))';
h = reshape(h,1,[]);
end
