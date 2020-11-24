function [ output ] = trimm_std( vector )

vals=prctile(abs(vector), 50);
output=std(vector(vector<vals & vector>-vals));

end

