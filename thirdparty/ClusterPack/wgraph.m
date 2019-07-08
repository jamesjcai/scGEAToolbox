% function dataname = wgraph(e,w,method,dataname)
%
% copyright (c) 1998-2011 by Alexander Strehl

function dataname = wgraph(e,w,method,dataname) 

if ~exist('method')
   method = 0;
end;

if ~exist('dataname'),
      dataname = ['' 'graph'];
end;
dataname = [dataname num2str(method)];

if ((method == 0)|(method == 1)),
  e=e-diag(diag(e));
end;
e = ints(e);
if ((method==1)|(method==3)),
  w = ints(w);
end;

while exist(dataname,'file')
   dataname = [dataname num2str(method)];
end;

fid = fopen(dataname,'w');

if (fid~=-1)

disp(['wgraph: writing ' dataname ]);

if (method == 0),
   fprintf(fid,'%d %d 1\n',[size(e,1) full(sum(sum(e>0)))/2]);
else
   if (method == 1),
      fprintf(fid,'%d %d 11\n',[size(e,1) full(sum(sum(e>0)))/2]);
   else
      validcolumns = find((sum(e,1)>0));
      if (length(validcolumns)~=size(e,2)),
         disp('wgraph: removing empty feature columns');
         e = e(:,validcolumns);
      end;
      fprintf(fid,'%d %d 1\n',[size(e,2) size(e,1)]);
   end;
end;

if ((method == 0)|(method == 1)),
   for i=1:size(e,1), 
      edges = find(e(i,:)>0);
      weights = e(i,edges);
      if method == 0,
         interlaced = zeros(1,2*length(edges));
         interlaced(1:2:2*length(edges)-1) = edges;
         interlaced(2:2:2*length(edges)) = weights;
      else
         interlaced = zeros(1,1+2*length(edges));
         interlaced(1) = w(i); 
         interlaced(2:2:2*length(edges)-1+1) = edges;
         interlaced(3:2:2*length(edges)+1) = weights;
      end;
      fprintf(fid,'%d ',interlaced); 
      fprintf(fid,'\n');
   end;
else
   disp(['wgraph: ' num2str(size(e,1)) ' vertices and ' num2str(size(e,2)) ' non-zero hyperedges']);
   for i=1:size(e,2), 
      edges = find(e(:,i)>0)';
      if method==2
         weight = sum(e(:,i));
      else
         weight = w(i);
      end;
      interlaced = [weight edges];
      fprintf(fid,'%d ',full(interlaced));
      fprintf(fid,'\n');
   end;
end;

fclose(fid);

else

disp(['wgraph: writing to ' dataname ' failed']);

end;

