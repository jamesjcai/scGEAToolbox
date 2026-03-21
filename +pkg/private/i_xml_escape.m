function s = i_xml_escape(s)
% I_XML_ESCAPE  Escape special characters for XML/HTML attribute values.
%   s = i_xml_escape(s) replaces &, <, >, " with their XML entity equivalents.
s = strrep(s, '&', '&amp;');
s = strrep(s, '<', '&lt;');
s = strrep(s, '>', '&gt;');
s = strrep(s, '"', '&quot;');
end
