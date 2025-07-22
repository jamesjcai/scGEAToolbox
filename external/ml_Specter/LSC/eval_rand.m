% Author: Modified package by Van Hoan Do
% Accuracy evaluation: returns the Rand Index given predicted and true
% labels
%
% Author: Frank Lin (frank@cs.cmu.edu)

function acc = eval_rand(truth, pred)

%i=ones(length(truth),1);

%acc=sum(sum((pred*i'==i*pred')==(truth*i'==i*truth')))/length(truth)^2;


acc = rand_index(truth, pred, 'adjusted');

end
