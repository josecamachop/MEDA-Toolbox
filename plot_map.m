
function fig_h = plot_map(map,label)

% Plot color map.
%
% plot_map(map) % minimum call
% plot_map(map,label) % complete call
%
%
% INPUTS:
%
% map: (MxM) matrix with values in the [0,1] interval. 
%
% label: (Mx1) name of the variables (numbers are used by default), eg.
%   num2str((1:M)')'
%
%
% OUTPUTS:
%
% fig_h: (1x1) figure handle.
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
%           Alejandro Perez Villegas (alextoni@gmail.com)
% last modification: 02/Feb/15.
%
% Copyright (C) 2014  University of Granada, Granada
% Copyright (C) 2014  Jose Camacho Paez
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

assert (nargin >= 1, 'Error: Missing arguments.');
assert (size(map,1) == size(map, 2), 'Dimension Error: map must be m-by-m.');
N = size(map,2);
if nargin < 2 || isempty(label)
    label=num2str((1:N)'); 
else
    if ischar(label), label = cellstr(label); end;
    assert (isequal(size(label), [N 1]), 'Dimension Error: label must be n-by-1.');
end

%% Main code

fig_h=figure;
map3 = [map map(:,end);map(end,:) map(end,end)];
sur_h=surface((1:N+1)'*ones(1,N+1),ones(N+1,1)*(1:N+1),map3);
if N < 100
    set(sur_h,'EdgeColor',[0.95 0.95 0.95]);
else
    set(sur_h,'EdgeColor','none');
end

% Set axis properties
axes_h = get(sur_h,'Parent');
set(axes_h,'Box','on');
set(axes_h,'XAxisLocation','top');
set(axes_h,'YDir','reverse');
set(axes_h,'XTick',(1:N)+0.5);
set(axes_h,'YTick',(1:N)+0.5);
set(axes_h,'XTickLabel',label);
set(axes_h,'YTickLabel',label);
if ~verLessThan('matlab', '8.3'),
    set(axes_h,'TicklabelInterpreter','None')
end

% Label font size
label_size = max(min(14,round(300/length(label))), 9);
set(axes_h, 'FontSize', label_size);

% Rotate X labels
rotateXLabels(axes_h,90);

% Resize axes position
pos = get(axes_h, 'Position');
set(axes_h,'Position',[pos(1) pos(2)/2 pos(3) pos(4)])

% Set colors
ind = [0:.2:0.79 0.8:0.04:1]';
set(fig_h,'Colormap',[[ind;ones(10,1)] [ind;flipud(ind)] [ones(10,1);flipud(ind)]])
caxis([-1 1]);
if find(map>0 & map<1)
    c_h=colorbar;
    set(c_h,'FontSize',14);
end

axis([1,N+1,1,N+1]);

    


        