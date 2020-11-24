function [ output ] = trimm_mean( vector )

% Di defalut rimuovo le code dello 0,25 % di valori

rimuovere=round(0.25/100*length(vector));

vector=sort(vector);
vector(1:rimuovere)=[];
vector(end-rimuovere:end)=[];


output=mean(vector);

end

