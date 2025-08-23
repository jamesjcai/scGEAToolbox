function h = i_toHash(txt)
% --- Helper: SHA-256 Hashing ---
    import java.security.*
    md = MessageDigest.getInstance('SHA-256');
    md.update(uint8(txt));
    d = typecast(md.digest, 'uint8');
    h = lower(dec2hex(d))';
    h = reshape(h,1,[]);
end