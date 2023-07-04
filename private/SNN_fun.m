function [A]=SNN_fun(d, outfile, k, distance)
% data input. row: object/sample; column: attribute.
% outfile: strength output
% k: nearest neighbor
% distance: methods for calculation distance (euclidean, correlation ...)

% source: https://github.com/BIOINSu/SNN-Cliq
% https://github.com/BIOINSu/SNN-Cliq/blob/master/Identification%20of%20cell%20types%20from%20single-cell%20transcriptomes%20using%20a%20novel%20clustering%20method%20.pdf

switch nargin
	case 1
		error('not enough arguments given: data, output file name.')
	case 2
		k=3;
		distance='euclidean';
	case 3
		distance='euclidean';
end

numSpl=size(d,1);
%%
IDX=knnsearch(d, d, 'K', k, 'Distance', distance);
strength=zeros(numSpl,numSpl);
for i=1:numSpl
for j=(i+1):numSpl
	nni=IDX(i,:);
	nnj=IDX(j,:);
	shared=intersect(nni,nnj);
    % the closeness depend on the rank of the shared knn in both list
	s=[0];	
	for l=1:length(shared)
		s=[s, k-0.5*(find(nni==shared(l))+find(nnj==shared(l)))];
	end
	% the score means how close the two points are.
	strength(i,j)=max(s);
	strength(j,i)=strength(i,j);
end
end
A=strength;
fn=fopen(outfile,'w');
for i=1:length(A)
	for j=1:i
		if(A(i,j)>0)
			fprintf(fn,'%d\t%d\t%f\r\n',i,j,A(i,j));
		end
	end
end
fclose(fn);

return