
function fig_h = plot_scatter(bdata,olabel,classes,axlabel,opt)

% Scatter plot.
%
% plot_scatter(bdata) % minimum call
% plot_scatter(bdata,olabel,classes,axlabel,opt) % complete call
%
%
% INPUTS:
%
% bdata: (Nx2) bidimensional data to plot. 
%
% olabel: {Nx1} name of the observations/variables (numbers are used by default)
%   use ' ' to avoid labels.
%
% classes: (Nx1) vector with the assignment of the observations/variables to classes, 
%   numbered from 1 onwards (1 class by default), eg. ones(N,1)
%
% axlabel: {2x1} variable/statistic plotted (nothing by default)
%
% opt: (1x1) options for data plotting.
%       0: filled marks (by default)
%       1: empty marks
%
%
% OUTPUTS:
%
% fig_h: (1x1) figure handle.
%
%
% coded by: José Camacho Páez (josecamacho@ugr.es).
% version: 2.1
% last modification: 03/Jul/14.
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
if nargin < 5, opt = 0; end;


%% Main code

if exist('classes')
    nc = max(classes);
    
    colors = ['b','g','r','c','m','y','k'];
    charac = ['o','*','s','d','v','^'];

    while length(colors)<nc,
        colors = [colors colors];
    end

    while length(charac)<nc,
        charac = [charac charac];
    end

else
    colors = char('b'*ones(1,s(1)));
    charac = char('o'*ones(1,s(1)));
    nc = 1;
    classes = ones(1,s(1));
end

if ~opt,
    colorsM = colors;
else
    colorsM = char('w'*ones(1,s(1)));
end
    
fig_h=figure;
hold on;
for i=1:nc,
    ind = find(classes==i);
    plot(bdata(ind,1),bdata(ind,2),strcat(colors(i),charac(i)),'MarkerFaceColor',colorsM(i));
end
if exist('axlabel')
    xlabel(axlabel{1},'FontSize',16);
    ylabel(axlabel{2},'FontSize',16);
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

    


        