function s=drawHierarchy(Sc,Tree,varargin)
% drawHierarchy Draw dendrogram for a hierarchical Tree
%
% Syntax
%__________________________________________________________________________
%
%   drawHierarchy(Sc,Tree)
%
%   drawHierarchy(Sc,Tree,Name,Value)
%
%
% Description
%__________________________________________________________________________
%
%   drawHierarchy(Sc,Tree) draws a dendrogram for a hierarchical cluster
%       Tree with initial partition Sc.
%
%   drawHierarchy(__,Name,Value) 
%
% 
% Input Arguments
%__________________________________________________________________________
%
%   Sc -- initial partition for the Tree
%
%   Tree -- Tree edges merging clusters in Sc
%
%
% Name-Value Pair Arguments
%__________________________________________________________________________
%
%   'Order' -- order of nodes to use for the visualization (should be
%       consistent with the Tree but this requirement is not checked)
%   
%   'Orientation' -- 'horizontal' (default), or 'vertical' 
%
%   'Axes' -- draw to 'Axes' instead of 'gca'
%
%   'LineWidth' -- width of dendrogram lines
%
%   'Color' -- color of dendrogram
%
%
% See also consensusPlot

% Version: 1.1.1
% Date: Thu  8 Mar 2018 15:34:45 CET
% Author: Lucas Jeub
% Email: ljeub@iu.edu

N=size(Sc,1);
parseArgs=inputParser();
addParameter(parseArgs,'Order',[]);
addParameter(parseArgs,'Orientation','vertical');
addParameter(parseArgs,'Axes',gca);
addParameter(parseArgs,'LineWidth',1);
addParameter(parseArgs,'Color',lines(1));
parse(parseArgs,varargin{:});

options=parseArgs.Results;

if ~ishold(options.Axes)
    cla;
end

if isempty(Tree)
    A=sparse(N,N);
    v_order=1;
else
% sort tree nodes
W=sparse(Tree(:,1),Tree(:,2),Tree(:,3),max(Sc),max(Sc));
A=sparse(Tree(:,1),Tree(:,2),1,max(Sc),max(Sc));
G=digraph(Tree(:,1),Tree(:,2));
v_order=dfsearch(G,1);
v_order=v_order(:)';
end
switch options.Orientation
    case 'horizontal'
        set(options.Axes,'ydir','reverse')
        draw_base=@(x,w,h) patch(options.Axes,[x-w/2,x+w/2,x],[1,1,1-h],options.Color,'LineStyle','none');
        draw_join=@(x1,x2,y1,y2) line(options.Axes,[x1,x1,x2],[y1,y2,y2],'color',options.Color,'LineWidth',options.LineWidth);
    case 'vertical'
        set(options.Axes,'xdir','reverse')
        draw_base=@(x,w,h) patch(options.Axes,[1,1,1-h],[x-w/2,x+w/2,x],options.Color,'LineStyle','none');
        draw_join=@(x1,x2,y1,y2) line(options.Axes,[y1,y2,y2],[x1,x1,x2],'color',options.Color,'LineWidth',options.LineWidth);
    otherwise
        error('unknown value for option orientation')
end

x=zeros(max(Sc),1);
y=ones(max(Sc),1);
if ~isempty(options.Order)
    x_base(options.Order)=1:N;
    s=options.Order;
else
    s=zeros(length(Sc),1);
end
pos=0;
for i=v_order(end:-1:1)
    ind=find(A(i,:));
    if isempty(ind)
        mm=find(Sc==i);
        if ~isempty(options.Order)
            x(i)=mean(x_base(mm));
        else
            x(i)=pos+1+length(mm)/2;
            s(pos+(1:length(mm)))=mm;
            pos=pos+length(mm);
        end
        parent=find(A(:,i));
        v=min(W(parent,A(parent,:)~=0));
        draw_base(x(i),length(mm),min((1-v),0.1));
        %patch([x(i)-length(mm)/2+0.5,x(i)+length(mm)/2-0.5,x(i)],[1,1,(y(i)+mean(A(A(:,i)~=0,i)))/2],lines(1));
    else
        x(i)=mean(x(ind));
        y(i)=min(W(i,ind));
        for j=ind(:)'
            draw_join(x(j),x(i),y(j),y(i));
            %line([x(j),x(j)],[y(j),y(i)]);
            %line([x(j),x(i)],[y(i),y(i)]);
        end
    end
end

switch options.Orientation
    case 'vertical'
        axis([0,1,0.5,N+0.5]);
        
    case 'horizontal'
        axis([0.5,N+0.5,0,1]);
end

end


