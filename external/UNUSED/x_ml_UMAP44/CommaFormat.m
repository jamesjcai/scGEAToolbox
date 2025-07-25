function [commaFormattedString] = CommaFormat(value)
%   AUTHORSHIP
%   Primary Developer:  Stephen Meehan <swmeehan@stanford.edu>
%   Math Lead & Secondary Developer: Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause

  % Split into integer part and fractional part.
  [integerPart, decimalPart]=strtok(num2str(value),'.'); 
  % Reverse the integer-part string.
  integerPart=integerPart(end:-1:1); 
  % Insert commas every third entry.
  integerPart=[sscanf(integerPart,'%c',[3,inf])' ... 
      repmat(',',ceil(length(integerPart)/3),1)]'; 
  integerPart=integerPart(:)'; 
  % Strip off any trailing commas.
  integerPart=deblank(integerPart(1:(end-1)));
  % Piece the integer part and fractional part back together again.
  commaFormattedString = [integerPart(end:-1:1) decimalPart];
  return; % CommaFormat
end