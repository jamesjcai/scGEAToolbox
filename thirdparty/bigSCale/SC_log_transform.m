function [ matrix ] = SC_log_transform( indices, total_data_norm , model )
% SC_log_transform ************************************************************************
% GIOVANNI IACONO, CNAG, 16/08/2017
% Transfomes a normalized expression matrix into different gene-normalized
% matrices
% INPUT
% indices: optional, restrict to "indices" genes. If empty uses all genes (Default).
% total_data_norm: normalized expression matrix
% model: see beloe for the descriptions
% OUTPUT 
% matrix: list of all significant markers




if isempty(indices)
    disp('I use all genes')
    indices=[1:length(total_data_norm(:,1))];
end

if (max(sum(total_data_norm))-min(sum(total_data_norm)))>1
    disp('Input data not normalized, you should pass normalized data!!');
    disp('Are you aware of this ?');
    pause(2);
end


switch model 
    
    case 1 % model=1. Log(x), then each row (gene) normalized between [-5:5]
    
        matrix=log2(total_data_norm(indices,:));
        dummy=unique(matrix);
        min_globale_detected=min(dummy(~isinf(dummy)));

        for k=1:length(indices)
            vettore=matrix(k,:);    
            vettore_inf=find(isinf(vettore));    
            vettore_num=find(~isinf(vettore));
            if length(unique(vettore(vettore_num)))>2
                x_out=shift_values(vettore(vettore_num),1,10);
            else
                x_out=zeros(size(vettore_num));
            end
            vettore(vettore_inf)=0;
            vettore(vettore_num)=x_out;
            matrix(k,:)=vettore;
        end

        nnz(isnan(matrix));
        matrix=matrix-5;

    case 2 % model=2. Log(x+1), then each row (gene) normalized between [-5:5]
     
        matrix=log2(total_data_norm(indices,:)+1);

	case 3 % model=3. sqrt(x), then each row (gene) normalized between [-5:5]
        
        matrix=sqrt(total_data_norm(indices,:));
            
	case 4 % model=4. Capped to 95% of expression. NO log nor sqrt, no normalization in [-5:5]
          matrix=total_data_norm;
          valori_estremi=prctile(matrix',95);
          valori_estremi_corr=prctile(matrix',99);
          valori_estremi(valori_estremi==0)=valori_estremi_corr(valori_estremi==0);
          for k=1:length(indices)
              matrix(k,matrix(k,:)>valori_estremi(k))=valori_estremi(k);
          end
              
end



if model>=2 && model<=4
    
	for k=1:length(indices)
        if length(unique(matrix(k,:)))>2
            x_out=shift_values(matrix(k,:),0,10);
        else
            x_out=zeros(size(matrix(k,:)));
        end
        matrix(k,:)=x_out;
    end
	matrix=matrix-5;
        
end

end
 


