
function fig_h = plot_Lscatter(bdata,olabel,classes,axlabel,opt,mult,maxv)

% Compressed scatter plot.
%
% plot_Lscatter(bdata) % minimum call
% plot_Lscatter(bdata,olabel,classes,axlabel),opt) % equivalent to plot_scatter
% plot_Lscatter(bdata,olabel,classes,axlabel),opt,mult,maxv) % complete call
%
%
% INPUTS:
%
% bdata: (Nx2) bidimensional data to plot. 
%
% olabel: {Nx1} name of the observations/variables (numbers are used by
%   default), use ' ' to avoid labels.
%
% classes: (Nx1) vector with the assignment of the observations/variables 
%   to classes, numbered from 1 onwards (1 class by default), eg. ones(N,1)
%
% axlabel: {2x1} variable/statistic plotted (nothing by default)
%
% opt: (1x1) options for data plotting.
%       0: filled marks, multiplicity information is not displayed
%       1: empty marks, multiplicity information is not displayed
%       2: 2D plot with the multiplicity info in the markers
%       3: 2D plot with the multiplicity info in the size of the markers
%           (by default)
%       4: 3D plot, with the multiplicity information in the Z axis
%       5: 2D hexagonal binning plot, with the multiplicity info in the 
%           size of the markers and class in the Z axis.
%
% mult: (Nx1) multiplicity of each row (1s by default)
%
% maxv: (1x3) thresholds for the different markers.
%       maxv(1): maximum threshold for marker 'x' for opt = 2 (20 by default)
%       maxv(1): maximum threshold for marker 'o' for opt = 2 (50 by default)
%       maxv(1): maximum threshold for marker 's' for opt = 2 (100 by default)
%
%
% OUTPUTS:
%
% fig_h: (1x1) figure handle.
%
%
% coded by: José Camacho Páez (josecamacho@ugr.es)
% last modification: 09/Apr/14.
%
% Copyright (C) 2014  José Camacho Páez
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
    
%

%% Parameters checking

if nargin < 1, error('Error in the number of arguments.'); end;
s = size(bdata);
if length(s) ~= 2 || s(2)~=2, error('Error in the dimension of the arguments.'); end;
if nargin < 2 || isempty(olabel)
    olabel=num2str((1:s(1))'); 
elseif ~isequal(olabel,' '),
    if ndims(olabel)==2 & find(size(olabel)==max(size(olabel)))==2, olabel = olabel'; end
    if size(olabel,1)~=s(1), error('Error in the dimension of the arguments.'); end;
end
if nargin < 3 || isempty(classes)
    classes = ones(s(1),1); 
else
    if ndims(classes)==2 & find(size(classes)==max(size(classes)))==2, classes = classes'; end
    if size(classes,1)~=s(1), error('Error in the dimension of the arguments.'); end;
end
if nargin < 4 ||isempty(axlabel)
    axlabel = {'Dim 1','Dim 2'}'; 
else
    if ndims(axlabel)==2 & find(size(axlabel)==max(size(axlabel)))==2, axlabel = axlabel'; end
    if size(axlabel,1)~=2, error('Error in the dimension of the arguments.'); end;
end
if nargin < 5, opt = 3; end;
if nargin < 6, mult = ones(s(1),1); end;
sm = size(mult);
if ndims(mult)==2 & find(size(mult)==max(size(mult)))==2, mult = mult'; end;
if s(1)~=sm(1) || sm(2)~=1 , error('Error in the dimension of the arguments.'); end;
if nargin < 7, maxv = [20 50 100]; end;
if ndims(maxv)==2 & find(size(maxv)==max(size(maxv)))==1, maxv = maxv'; end;
if size(maxv,2)~=3 || size(maxv,1)~=1, error('Error in the dimension of the arguments.'); end;


%% Main code

chr = ['.','x','o','s','d'];
maxv = [0 1 maxv Inf];

if exist('classes')
    m = max(classes);
    
    colors = ['b','g','r','c','m','y','k'];
    charac = ['o','*','s','d','v','^'];

    while length(colors)<m,
        colors = [colors colors];
    end

    while length(charac)<m,
        charac = [charac charac];
    end

else
    colors = char('b'*ones(1,s(1)));
    charac = char('o'*ones(1,s(1)));
    m = 1;
    classes = ones(1,s(1));
end
    
fig_h=figure;
hold on;

for i=1:m,
    ind=find(classes==i);
    switch opt,
        case 2,
             for j=1:length(maxv)-1,
                ind2 = find(mult(ind)<=maxv(j+1));
                ind3 = find(mult(ind(ind2))>maxv(j));
                ind2 = ind(ind2(ind3));
                plot(bdata(ind2,1),bdata(ind2,2),strcat(colors(i),chr(j)),'MarkerFaceColor',colors(i));
            end
            %if i==1, legend('n=1',['1<n\leq' sprintf('%d',maxv(3))],[sprintf('%d',maxv(3)) '<n\leq' sprintf('%d',maxv(4))],[sprintf('%d',maxv(4)) '<n\leq' sprintf('%d',maxv(5))],sprintf('%d<n',maxv(5)),'Location','Best'); end   
        case 3,
             for j=1:length(maxv)-1,
                ind2 = find(mult(ind)<=maxv(j+1));
                ind3 = find(mult(ind(ind2))>maxv(j));
                ind2 = ind(ind2(ind3));
                plot(bdata(ind2,1),bdata(ind2,2),strcat(colors(i),charac(i)),'MarkerFaceColor',colors(i),'MarkerSize',round(j*2.5));
            end
            %if i==1, legend('n=1',['1<n\leq' sprintf('%d',maxv(3))],[sprintf('%d',maxv(3)) '<n\leq' sprintf('%d',maxv(4))],[sprintf('%d',maxv(4)) '<n\leq' sprintf('%d',maxv(5))],sprintf('%d<n',maxv(5)),'Location','Best'); end   
        case 4,
            for j=1:length(ind),
                plot3(bdata(ind([j j]),1),bdata(ind([j j]),2),[0 mult(ind(j))],'-','color',colors(i),'MarkerFaceColor',colors(i));
                plot3(bdata(ind(j),1),bdata(ind(j),2),mult(ind(j)),strcat(colors(i),charac(i)),'MarkerFaceColor',colors(i));
            end
        case 5,
             for j=1:length(maxv)-1,
                ind2 = find(mult(ind)<=maxv(j+1));
                ind3 = find(mult(ind(ind2))>maxv(j));
                ind2 = ind(ind2(ind3));
                plot3(bdata(ind2,1),bdata(ind2,2),i*ones(length(ind2),1),strcat(colors(i),charac(i)),'MarkerFaceColor',colors(i),'MarkerSize',round(j*2.5));
            end
        otherwise
            if ~opt,
                colorsM = colors;
            else
                colorsM = char('w'*ones(1,s(1)));
            end
            plot(bdata(ind,1),bdata(ind,2),strcat(colors(i),charac(i)),'MarkerFaceColor',colorsM(i));
    end
    if exist('axlabel')
        xlabel(axlabel{1},'FontSize',16);
        ylabel(axlabel{2},'FontSize',16);
    end

    labels{i}=sprintf('class %d',i);
end
        
ax=axis;
plot([0 0],ax(3:4),'k--');
plot(ax(1:2),[0 0],'k--');
axis(ax);
axes_h=get(fig_h,'Children');
axes_h=axes_h(1);
set(axes_h,'FontSize',14);
if ~isequal(olabel,' '),  
    deltax = (ax(2) - ax(1))/50;
    deltay = (ax(4) - ax(3))/50;
    text(bdata(:,1)+deltax,bdata(:,2)+deltay,olabel);
end
box on

    


        