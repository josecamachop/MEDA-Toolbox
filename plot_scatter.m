
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
% olabel: {Nx1} name of the observations/variables
%   Allowed cell array of strings, eg. {'first', 'second', 'third', ...}
%   use [] to set the default, empty labels.
%
% classes: (Nx1) vector with the assignment of the observations/variables to classes,
%   Allowed numerical classes, eg. [1 1 2 2 2 3], 
%   and cell array of strings, eg. {'blue','red','red','green','blue'}.
%   use [] to set the default, a single class.
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
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
%           Alejandro Perez Villegas (alextoni@gmail.com)
% last modification: 25/Jan/15.
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
N = size(bdata, 1);
if nargin < 2 || isempty(olabel)
    olabel = repmat({''}, N, 1);
end
if nargin < 3 || isempty(classes)
    classes = ones(N, 1);
end
if nargin < 4 || isempty(axlabel)
    axlabel = {'Dim 1';'Dim 2'};
end
if nargin < 5, opt = 0; end;

% Convert char arrays to cell
if ischar(olabel),  olabel = cellstr(olabel); end;
if ischar(classes), classes = cellstr(classes); end;
if ischar(axlabel), axlabel = cellstr(axlabel); end;

% Convert row arrays to column arrays
if size(olabel,1)  == 1, olabel = olabel'; end;
if size(classes,1) == 1, classes = classes'; end;
if size(axlabel,1) == 1, axlabel = axlabel'; end;

% Validate dimensions of input data
assert (size(bdata,2) == 2, 'Dimension Error: bdata must be n-by-2.')
assert (isequal(size(olabel), [N 1]), 'Dimension Error: label must be n-by-1.');
assert (isequal(size(classes), [N 1]), 'Dimension Error: classes must be n-by-1.')
assert (isequal(size(axlabel), [2 1]), 'Dimension Error: axlabel must be 2-by-1.')

%% Main code
fig_h = figure;
hold on;

% Plot points
unique_classes = unique(classes);
for i=1:length(unique_classes)
    ind = find(classes == unique_classes(i));
    a= scatter(bdata(ind,1), bdata(ind,2), [], 'o');
end

% Fill marks
if opt == 0
    for i=1:length(a)
        color = get(a(i), 'Color');
        set(a(i), 'MarkerFaceColor',color);
    end
end

% Plot labels
text(bdata(:,1), bdata(:,2), olabel(:,1), 'VerticalAlignment','bottom', 'HorizontalAlignment','left');

% Set axis labels
xlabel(axlabel(1), 'FontSize', 16);
ylabel(axlabel(2), 'FontSize', 16);

% Plot origin lines
ax = axis;
plot([0 0], ax(3:4), 'k--', 'HandleVisibility', 'off');
plot(ax(1:2), [0 0], 'k--', 'HandleVisibility', 'off');
axis(ax)

legend off
box on

