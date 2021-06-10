function i_markergenesviolin(X,genelist,c,markerlist,numfig,N,sptltxt)

if nargin<7, sptltxt=''; end
if nargin<6, N=4; end
n=length(markerlist);
if nargin<5 || isempty(numfig)
    numfig=ceil(n/N);
end
if ~ismember(N,[4 9 16]) || numfig<1 || numfig>10
    error('error');
end

for kkk=1:numfig
    figurec;
    for kk=1:min([N,n-(kkk-1)*N])
        if N==16
            subplot(4,4,kk);
        elseif N==9
            subplot(3,3,kk)
        elseif N==4
            subplot(2,2,kk);
        end
        currentg=markerlist(kk+N*(kkk-1));
        idx=genelist==currentg;
        x=log2(X(idx,:)+1);
        pkg.i_violinplot(x,c);
        xtickangle(-45);
        title(currentg);
    end
    if ~isempty(sptltxt)
        suptitle(sptltxt);
    end
end
end


% FIGUREC - create a figure window in a non-overlapping (cascading)
%           location
% 
% USAGE:
% 
% figurec
% figurec(...)
% h=figurec
% h=figurec(...)
% 
% FIGUREC acts just like the Matlab FIGURE command, with all arguments
% passed through, except that the new figure is created a little to the
% right and down from the highest numbered figure currently existing, so
% that they won't overlap. If moving the location would push the figure too
% close to the edge of the screen, then the new figure is created in the
% default location as usual. (Subsequent figures will again be cascaded.)
% 
% EXAMPLE:
% 
% close all;for n=1:20;figurec('color',rand(1,3));plot(0,0);title('Sample');end
function varargout=figurec(varargin)
f = findobj(0,'type','figure'); % list of existing figures
ss=get(0,'ScreenSize'); % pixel size of entire screen
h=figure(varargin{:}); % create figure using pass-through arguments
hp = get(h,'pos'); % size of new figure when created
if ~isempty(f)
    f = f(1);
    u=get(f,'units');
    set(f,'units','pixels')
    p=get(f,'pos');
    set(f,'units',u)
    % if moving won't push too far, move; else leave in default location
    if p(1)+50+hp(3) <= ss(3) && p(2) >= 5
        u=get(h,'units');
        ss = get(0,'screensize');
        set(h,'units','pixels')
        set(h,'pos',[p(1)+50 p(2)-50 hp(3:4)]);
        set(h,'units',u)
    end
end
if nargout > 0
    varargout{1}=h;
end
end

