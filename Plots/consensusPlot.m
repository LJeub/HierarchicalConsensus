function [ax_C,ax_H,order]=consensusPlot(C,Sc,Tree,varargin)
% consensusPlot Plot hierarchical consensus clustering results
%
% Syntax
%__________________________________________________________________________
%
%   consensusPlot(C,Sc,Tree)
%
%   consensusPlot(__,Name,Value)
%
%   [ax_C,ax_H]=consensusPlot(__)
%
%
% Description
%__________________________________________________________________________
%
%   consensusPlot(C,Sc,Tree) draws coclassification matrix and consensus
%       hierarchy
%
%   consensusPlot(__,Name,Value) customizes the plot based on 'Name-Value'
%       arguments.
%
%   [ax_C,ax_H]=consensusPlot(__) returns the axes used for plotting
%       (useful for later customization).
%
%
% Input Arguments
%__________________________________________________________________________
%
%   C -- Coclassification Matrix.
%
%   Sc -- Initial partition for Tree.
%
%   Tree -- Hierarchical cluster Tree for merging clusters in 'Sc'
%
%
% Name-Value Pair Arguments
%__________________________________________________________________________
%
%   'GroundTruth' -- Optionally provide ground-truth partition(s) (shown as
%       background colors under hierarchy plot).
%
%   'Labels' -- Node labels (shown on y-axis)
%
%   'LabelFont' -- Font size for node labels (default: 6)
%
%   'ImageRescale' -- Rescale images (useful for printing to pdf to avoid
%       blurring in pdf viewers that interpolate) (default: 1)
%
% Output Arguments
%__________________________________________________________________________
%
%   ax_C -- Axes used to draw coclassification matrix
%
%   ax_H -- Axes used to draw consensus hierarchy
%
% See also hierachicalConsensus, drawHierarchy

% Version: 1.0
% Date: Thu Oct  5 17:09:19 EDT 2017
% Author: Lucas Jeub
% Email: ljeub@iu.edu


parseArgs=inputParser();
addParameter(parseArgs,'GroundTruth',[]);
addParameter(parseArgs,'Labels',{});
addParameter(parseArgs,'LabelFont',6);
addParameter(parseArgs,'ImageRescale',1);
parse(parseArgs,varargin{:});
Sgtruth=parseArgs.Results.GroundTruth;
labels=parseArgs.Results.Labels;
ylabelfont=parseArgs.Results.LabelFont;
imrescale=parseArgs.Results.ImageRescale;
clf()
ax_C_int=axes('position',[0.1,0.1,0.7,0.8]);
ax_H_int=axes('position',[0.8,0.1,0.1,0.8]);

N=length(C);
ylims=[0.5+0.5*1/imrescale,N+0.5-0.5*1/imrescale];
s_int=treeSort(C,Sc,Tree);
imagesc(ax_C_int,ylims,ylims, imresize(C(s_int,s_int),imrescale,'nearest'));
if ~isempty(labels)
    set(ax_C_int,'ytick',1:N,'yticklabel',labels(s_int),...
        'xtick',1:N,'xticklabel',labels(s_int),'xticklabelrotation',90,...
        'ticklabelinterpreter','None','fontsize',ylabelfont)
end

if ~isempty(Sgtruth)
w=1/size(Sgtruth,2);
xlims=[w/2*1/imrescale,1-w/2*1/imrescale];

imdata=zeros(N,size(Sgtruth,2),3);
for j=1:size(Sgtruth,2)
    if max(Sgtruth(:,j))<=6
        gtruthcolors=lines(max(Sgtruth(:)));
    else
        gtruthcolors=parula(max(Sgtruth(:)));
        gtruthcolors=gtruthcolors(randperm(size(gtruthcolors,1)),:);
    end
    for i=1:N
        imdata(i,j,:)=gtruthcolors(Sgtruth(s_int(i),j),:);
    end
end
image(ax_H_int,xlims,ylims,imresize(imdata,imrescale,'nearest'),'AlphaData',0.3);
set(ax_H_int,'nextplot','add')
end
drawHierarchy(Sc,Tree,'axes',ax_H_int,'order',s_int,'orientation','vertical');
set(ax_H_int,'ytick',[],'ydir','reverse');
if nargout>0
    ax_C=ax_C_int;
    ax_H=ax_H_int;
    order=s_int;
end
end

