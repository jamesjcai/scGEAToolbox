%
FACTOR=max_factor;
decrease=0;
all_cells=uint32(1:length(cells));
cells_in_use=all_cells;

picked=randi(numel(distances),20000000,1);
MAX_DIST=prctile(distances(picked), [5 10 20 30 40 50]);
%MAX_DIST=prctile(reshape(distances,numel(distances),1),[5 10 20 30 40 50]);
disp('Calcolated percentile')


% ordino le cellule per distanza
for k=1:length(distances(:,1))
    [c ix]=sort(distances(k,:));
    distances(k,:)=c;
    cells(k,:)=cells(k,ix);
end
clear c;
disp('Ordered cells ')
blocks=uint32([]);

STEP_DIST=1;





max_dist=MAX_DIST(STEP_DIST)

conta=1;

paired_cells_old=uint32([0]);
round=0;
while 1
    round=round+1;
    while 1
        if FACTOR>1   
            available=find ( sum(~isinf(distances(:,1:FACTOR)')) == FACTOR);
        else
            available=find ( (~isinf(distances(:,1:FACTOR)')) == FACTOR);
        end

        if any(available)
            [val ix]=min(distances(available,FACTOR));
            ix=available(ix);
            if (isinf(val) | val>max_dist)
                break;
            end
            blocks=[ blocks; cells_in_use(ix) cells(ix,1:FACTOR) nan(1,decrease)];
            
            converted_cells=find(ismember(cells_in_use,blocks(end,:)));
            distances(converted_cells,:)=Inf; % metto a Inf le righe delle cellule usate
            distances(ismember(cells(:,1:FACTOR),blocks(end,:)))=Inf; % metto a Inf le celle giuste nella colonna factor
            
            conta=conta+1;
        else 
            break;
        end
    end
    
    paired_cells=nonzeros(reshape(blocks,numel(blocks),1));
    unpaired_cells=setdiff(all_cells,paired_cells);
    sprintf('Round %g, covered %g cells, block size =%g,max_dist=%.1f (STEP_DIST=%g)',round,nnz(~isnan(paired_cells)),FACTOR+1,max_dist,STEP_DIST)

    % VELOCIZZABILE FACENDO IS MEMBER SOLO SULLE CELLULE CHE NON SONO STATE
    % GIA RIMOSSE
    converted_unpaired=find(ismember(cells_in_use,unpaired_cells));
    
    if round<5
        save('distances_temp.mat','distances','-v7.3');
        clear distances
        cells=cells(converted_unpaired,:);
        load('distances_temp.mat');
    else
        cells=cells(converted_unpaired,:);
    end
    
    distances=distances(converted_unpaired,:);
    cells_in_use=cells_in_use(converted_unpaired);
    distances(ismember(cells,setdiff(paired_cells,paired_cells_old)))=Inf;
    
    % ordino le cellule per distanza
    for k=1:length(distances(:,1))
        [c ix]=sort(distances(k,:));
        distances(k,:)=c;
        cells(k,:)=cells(k,ix);
    end
    clear c;

    if numel(paired_cells)==numel(paired_cells_old)
        disp('Converged');
        
        FACTOR=FACTOR-1;
        decrease=decrease+1;
        sprintf('Reducing factor to %g',FACTOR)
        
        if ( FACTOR==0)% | STEP_DIST>adm_max(FACTOR) )
            
            if (STEP_DIST==length(MAX_DIST))
                disp('Nothing more I can do');
                return;
            end
            
            FACTOR=max_factor;
            decrease=0;
            STEP_DIST=STEP_DIST+1;
            max_dist=MAX_DIST(STEP_DIST);
            
        else
            max_dist=MAX_DIST(STEP_DIST);
        end

    else
        paired_cells_old = paired_cells;
    end
end











%     if numel(paired_cells)==numel(paired_cells_old)
%         disp('Arrivati a convergenza');
%         
%         STEP_DIST=STEP_DIST+1;
%        
%         if ( STEP_DIST>adm_max(FACTOR) | STEP_DIST>length(MAX_DIST) )
%             
%             if FACTOR==1    
%                 disp('Siamo alla frutta');
%                 return;
%             end
%             
%         STEP_DIST=1;
%         max_dist=MAX_DIST(STEP_DIST);
%         FACTOR=FACTOR-1;
%         decrease=decrease+1;
%         sprintf('Riduco factor a %g',FACTOR)
%         
%         else
%             max_dist=MAX_DIST(STEP_DIST);
%         end
% 
%     else
%         paired_cells_old = paired_cells;
%     end



% function [cells distances cells_paired cells_used clones final Closests Distances]=SC_bigdata_collapse( distances, cells )
% 
% 
% % ordino le cellule per distanza
% [c ix]=sort(distances,2);
% for k=1:length(cells(:,1))
%         cells(k,:)=cells(k,ix(k,:));
% end
% distances=c;
% clear c;
% 
% 
%     
% for h=1:5
%     h
%     [cells distances]=collapse_fast( distances, cells );
%     nnz(unique(cells(:,1)))
% end
% 
% 
% 
% cells_paired=zeros(length(cells),1);
% cells_used=zeros(length(cells),1);
% clones=zeros(length(cells),1);
% 
% for h=1:length(cells)
%     if   nnz(cells(:,1)==cells(h,1))==1
%     cells_paired(h)=1;
%     cells_used(cells(h,1))=1;
%     end
%     
%     %dummy=cells(setdiff(1:length(cells),h),1); % prima colonna meno la riga in questione    
%     %clones(h)=nnz(unique(dummy(ismember(dummy,cells(h,:)))));
% end
% 
% 
% 
% cells_unpaired=find(cells_paired==0);
% cells_sub=cells(cells_unpaired,:);
% distances_sub=distances(cells_unpaired,:);
% 
% % rimuovo le cellule che sono state usate
% used=find(ismember(cells_sub,find(cells_used)));
% cells_sub(used)=0;
% distances_sub(used)=Inf;
% 
% final=zeros(length(cells),1);
% [ Closests Distances sub_final]=collapse_slow( distances_sub, cells_sub );
% final(cells_unpaired)=sub_final;
% 
% 
% 
% 
% 
% end
% 
% 
% 
% 
% 
% 
% 
% 
% function [ Closests Distances final_assignment]=collapse_slow( distances, cells )
% 
% 
%     dims=size(cells);
% 
%     Closests=zeros(length(cells),1);
%     Distances=zeros(length(cells),1);
%     Counts=ones(length(cells),1);
%     all_cells=unique(cells);
% 
%     for k=1:length(all_cells)
%         
%         pos=find(cells==all_cells(k));
%         dummy=distances(pos);
%         [c ix]=sort(dummy);
%         [closest_r closest_c]= ind2sub( dims, pos(ix(1)));
%         Closests(closest_r,Counts(closest_r))=all_cells(k);
%         Distances(closest_r,Counts(closest_r))=c(1);
%         Counts(closest_r)=Counts(closest_r)+1;
%     end
% 
%     final_assignment=zeros(length(cells(:,1)),1);
% 
%     for k=1:length(cells)
%         
%         [c ix]=min(Distances(k,1:Counts(k)-1));
%         if ~isempty(c)
%            final_assignment(k)=Closests(k,ix(1));
%         end
%     end
% 
% end
% 
% 
% function [cells distances]=collapse_fast( distances, cells )
% 
% size_cells=length(cells(1,:))
% 
%     for k=1:length(cells)
%         
%         conta=1;
%         add=0;
%         while conta<=size_cells
%             %disp('**************************************************')
%             %k
%             %conta  
%             %cells(k,conta)
%             %distances(k,conta)
%             pos=find(cells(:,1)==cells(k,conta));
%             %min(distances(pos,1))
%             %input('')
%             if (length(pos)+add)>1
%                     if distances(k,conta)<=min(distances(pos,1))
%                         dummy=cells(k,1);
%                         cells(k,1)=cells(k,conta);
%                         cells(k,conta)=dummy; 
%                         dummy=distances(k,1);
%                         distances(k,1)=distances(k,conta);
%                         distances(k,conta)=dummy;
%                         break;
%                     else
%                         conta=conta+1;
%                         add=1;
%                     end
%             else
%                 dummy=cells(k,1);
%                 cells(k,1)=cells(k,conta);
%                 cells(k,conta)=dummy; 
%                 dummy=distances(k,1);
%                 distances(k,1)=distances(k,conta);
%                 distances(k,conta)=dummy;
%                 break;
%             end % if length(pos)
%         end %while
%         if conta==size_cells
%             sprintf('k=%g)Raggiunto conta %g',k,conta)
%         end
%     end
%     
% end


