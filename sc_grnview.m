function [hs]=sc_grnview(A,g,alpha)

if nargin<3, alpha=0.75; end
n=size(A,1);
if n>200, error('Smaller network is expected'); end
if nargin<2, g=string(1:n); end

% Add the UI components
hs = addcomponents;
[G,p]=drawnetwork(A,g,alpha);   
% Make figure visible after adding components
hs.fig.Visible = 'on';
   
   function hs = addcomponents
       % Add components, save handles in a struct
       hs.fig = figure('Visible','off',...
                  'Tag','fig',...
                  'SizeChangedFcn',@resizeui);
              
       hs.btn = uicontrol(hs.fig,'String',...
                  'Layout',...
                  'Callback',@ChageLayout,...
                  'Tag','button');
              
       hs.btn2 = uicontrol(hs.fig,'String',...
                  'Weight',...
                  'Callback',@ChageWeight,...
                  'Tag','button');   
              
       hs.btn3 = uicontrol(hs.fig,'String',...
                  'Cutoff',...
                  'Callback',@ChageCutoff,...
                  'Tag','button');   
              
       hs.btn4 = uicontrol(hs.fig,'String',...
                  'Directed',...
                  'Callback',@ChageDirected,...
                  'Tag','button');   
       hs.btn5 = uicontrol(hs.fig,'String',...
                  'Font Size',...
                  'Callback',@ChageFontSize,...
                  'Tag','button');   
              
       hs.ax = axes('Parent',hs.fig,...
                  'Units','pixels',...
                  'Tag','ax');
   end

   function ChageLayout(hObject,event)
       %plot(hs.ax,theta,y);
       a=["layered","subspace","force","circle","force3","subspace3"];       
       i=randi(length(a));       
       p.layout(a(i));
       if i>4
           view(3);
       else
           view(2);
       end
   end

   function ChageWeight(hObject,event)
       if length(unique(p.LineWidth))>1
        p.LineWidth = p.LineWidth./p.LineWidth;
       else
        a=3:10;
        b=a(randi(length(a),1));
        G.Edges.LWidths = abs(b*G.Edges.Weight/max(G.Edges.Weight));
        p.LineWidth = G.Edges.LWidths;
       end
   end

   function ChageCutoff(hObject,event)
        list = {'0.00 (show all edges)',...
            '0.30','0.35','0.40','0.45',...
            '0.50','0.55','0.60',...
            '0.65','0.70','0.75','0.80','0.85',...
            '0.90','0.95 (show 5% of edges)'};
        [indx,tf] = listdlg('ListString',list,...
            'SelectionMode','single','ListSize',[160,230]);
        if tf
            if indx==1
                cutoff=0;
            elseif indx==length(list)
                cutoff=0.95;
            else
                cutoff=str2double(list(indx));
            end
            [G,p]=drawnetwork(A,g,cutoff);
        end
   end

   function ChageDirected(hObject,event)
        if issymmetric(G.adjacency)            
            [G,p]=drawnetwork(A,g,0.65);
        else
            [G,p]=drawnetwork(0.5*(A+A.'),g,0.65);
        end
   end

   function ChageFontSize(hObject,event)
       if p.NodeFontSize>=20
           p.NodeFontSize=7;
       else
           p.NodeFontSize=p.NodeFontSize+1;
       end
   end

   function resizeui(hObject,event)
           
       % Get figure width and height
       figwidth = hs.fig.Position(3);
       figheight = hs.fig.Position(4);
       
       % Set button position
       bheight = 30; 
       bwidth = 70;
       bbottomedge = figheight - bheight - 50;
       bleftedge = 10;
       hs.btn.Position = [bleftedge bbottomedge bwidth bheight];
       hs.btn2.Position = [bleftedge bbottomedge-40 bwidth bheight];
       hs.btn3.Position = [bleftedge bbottomedge-80 bwidth bheight];
       hs.btn4.Position = [bleftedge bbottomedge-120 bwidth bheight];
       hs.btn5.Position = [bleftedge bbottomedge-160 bwidth bheight];
       
       % Set axes position
       axheight = .85*figheight;
       axbottomedge = max(0,figheight - axheight - 30);
       axleftedge = bleftedge + bwidth + 10;
       axwidth = max(0,figwidth - axleftedge - 30);
       hs.ax.Position = [axleftedge axbottomedge axwidth axheight];
   end

    function [G,p]=drawnetwork(A,g,alpha)
        B=e_transf(A,alpha);
        if issymmetric(B)
            G=graph(B,g,'omitselfloops');    
        else
            G=digraph(B,g,'omitselfloops');
        end
        p=plot(hs.ax,G);
        layout(p,'force');
        if isa(G,'digraph')
            G.Nodes.NodeColors = outdegree(G)-indegree(G);
        else
            G.Nodes.NodeColors = degree(G);            
        end
        p.NodeCData = G.Nodes.NodeColors;
        n=size(G.Edges,1);
        cc=repmat([0 0.4470 0.7410],n,1);
        cc(G.Edges.Weight<0,:)=repmat([0.8500, 0.3250, 0.0980],...
               sum(G.Edges.Weight<0),1);
        p.EdgeColor=cc;
        p.NodeFontSize=2*p.NodeFontSize;
        %view(3)            
        title('Single-cell Gene Regulatory Netowrk');
    end
end


function A=e_transf(A,q)
% A - adjacency matrix
if nargin<2, q=0.95; end
dim=size(A);
if numel(dim)==2
    a=max(abs(A(:)));
    if a>0
        A=A./a;
        A=A.*(abs(A)>quantile(abs(A(:)),q));        
    end
elseif numel(dim)==3
    for k=1:dim(3)
        A(:,:,k)=e_transf(A(:,:,k),q);
    end
end
end