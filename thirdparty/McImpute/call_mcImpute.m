% function [Xrec,normalizedXrec,geneGilteredXrec,genelist] = call_mcImpute(X,genelist)
% calling fucntion for mcImpute

% MANDATORY ARGUMENT-> data: scRNA seq data matrix (cell x gene)
% 
%OPTIONAL ARGUMENT-> dataname: string having dataname, to save library size
%                    of data with unique name
%OPTIONAL ARGUMENT-> gene_names, gene_ids
%OPTIONAL ARGUMENT-> pro_dir: path to directory where processed matrix
%                    should be saved, if not specified a new directory would be created in
%                    current folder and processed data would be saved their

%% variables initialization
% dataname='NEWDATA'; gene_names=[]; gene_ids=[]; pro_dir='./processed_data/';
% 
% if ~isempty(varargin)
%     for j = 1:length(varargin)
%        
%         if strcmp(varargin{j}, 'dataname')
%             dataname = varargin{j+1};
%         end
%     
%         if strcmp(varargin{j}, 'gene_names')
%             gene_names = varargin{j+1};
%         end
%         
%         if strcmp(varargin{j}, 'gene_ids')
%             gene_ids = varargin{j+1};
%         end
%         
%         if strcmp(varargin{j}, 'pro_dir')
%             pro_dir = varargin{j+1};
%         end
%         
%     end
% end

%{ 
elseif (nargin<2)
    dataname='NEWDATA'; gene_names=[]; gene_ids=[]; pro=1;
    warning('Dataname not provided..library size saved with prefix NEWDATA..')
elseif(nargin<3)
    gene_names=[]; gene_ids=[]; pro=1;
elseif(nargin<4)
    gene_ids=1:length(gene_names); pro=1;
else
    pro=1;
end
%}



% if ~exist('Temp','dir'), mkdir Temp; end
% if ~exist(pro_dir,'dir'), mkdir processedData; pro_dir='processedData/'; end


%% Process data
% pro=1;
% if (pro)
% disp('Processing data.....')
% [processed_data,gene_names] = process(data,dataname,gene_names, gene_ids); %processed_data=data
% % save([pro_dir '' dataname '_processed.mat'],'processed_data')
% end

[X,genelist]=sc_selectg(X,genelist,2,3);
libsize=sum(X)';
X=sc_norm(X);
X=log2(1+X);

% a2=log2(1+sc_norm(sc_selectg(X',genelistx,2,3)')');
processed_data=X';

%% Run mcimpute
IDX = find(processed_data);
M = opRestriction(numel(processed_data),IDX);
y = M(processed_data(:),1);
 
[Xrec] = IST_MC(y,M,size(processed_data),0); %mask changed and lansvd changed, with NN constraint
%Xrec = IST_eMC(y,M,size(processed_data),11); 
normalizedXrec=2.^(Xrec)-1;
normed_data=2.^Xrec-1;
raw_data=(normed_data.*libsize)./ median(libsize);
geneGilteredXrec=round(raw_data);

