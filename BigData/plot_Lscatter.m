
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
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
%           Alejandro Perez Villegas (alextoni@gmail.com)
% last modification: 06/Feb/15.
%
% Copyright (C) 2016  University of Granada, Granada
% Copyright (C) 2016  Jose Camacho Paez
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
% Set default values
assert(nargin >= 1, 'Error: Missing arguments.');
N = size(bdata, 1);
if nargin < 2 || isempty(olabel),  olabel  = repmat({''},N,1);  end;
if nargin < 3 || isempty(classes), classes = ones(N, 1);        end;
if nargin < 4 || isempty(axlabel), axlabel = {'Dim 1';'Dim 2'}; end;
if nargin < 5 || isempty(opt),     opt     = 3;                 end;
if nargin < 6 || isempty(mult),    mult    = ones(N,1);         end;
if nargin < 7 || isempty(maxv),    maxv    = [20 50 100];       end;

% Convert char arrays to cell
if ischar(olabel),  olabel  = cellstr(olabel);  end;
if ischar(classes), classes = cellstr(classes); end;
if ischar(axlabel), axlabel = cellstr(axlabel); end;

% Convert row arrays to column arrays
if size(olabel,1)  == 1, olabel  = olabel';  end;
if size(classes,1) == 1, classes = classes'; end;
if size(axlabel,1) == 1, axlabel = axlabel'; end;
if size(mult,1) == 1, mult = mult'; end;
if size(maxv,2) == 1, maxv = maxv'; end;

% Validate dimensions of input data
assert(size(bdata,2) == 2, 'Dimension Error: bdata must be n-by-2.');
assert(isequal(size(olabel),  [N 1]), 'Dimension Error: label must be n-by-1.');
assert(isequal(size(classes), [N 1]), 'Dimension Error: classes must be n-by-1.');
assert(isequal(size(axlabel), [2 1]), 'Dimension Error: axlabel must be 2-by-1.');
assert(isequal(size(mult), [N 1]), 'Dimension Error: mult must be n-by-1');
assert(isequal(size(maxv), [1 3]), 'Dimension Error: maxv must be 1-by-3');


%% Main code
% Preprocess classes to force them start with 1, 2...n,
unique_classes = unique(classes);
if iscell(classes)
    normal_classes = arrayfun(@(x) find(strcmp(unique_classes, x), 1), classes);
else
    normal_classes = arrayfun(@(x) find(unique_classes == x, 1), classes);
end

% Define mult bins, markers, colors and sizes 
bins = [0 1 maxv Inf];
markers = ['^','v','d','o','s'];

color_list = hsv(length(unique_classes));
colors = color_list(normal_classes, :);

sizes = zeros(size(mult));
for i=1:length(bins)-1
    sizes (mult>bins(i) & mult<=bins(i+1)) = round(2.5 * i^2 * pi);
end
 
% Plot points and labels
fig_h = figure;
hold on;
switch opt
    case 0,  % 2D plot, No multiplicity info, filled marks
        for i=1:length(unique_classes)
            ind = classes == unique_classes(i);
            scatter(bdata(ind,1), bdata(ind,2), [], colors(ind,:),'filled','DisplayName',num2str(unique_classes(i)));
        end
        text(bdata(:,1), bdata(:,2), olabel(:,1), 'VerticalAlignment','bottom','HorizontalAlignment','left');

    case 1,  % 2D plot, No multiplicity info, empty marks
        for i=1:length(unique_classes)
            ind = classes == unique_classes(i);
            scatter(bdata(ind,1), bdata(ind,2), [], colors(ind,:),'DisplayName',num2str(unique_classes(i)));
        end
        text(bdata(:,1), bdata(:,2), olabel(:,1), 'VerticalAlignment','bottom','HorizontalAlignment','left');
    
    case 2,  % 2D plot, Multiplicity in markers
        for i=1:length(unique_classes)
            for j=1:length(bins)-1
                ind = classes==unique_classes(i) & mult<=bins(j+1) & mult>bins(j);
                disp_name = strcat(num2str(unique_classes(i)), ' (mult: > ', num2str(bins(j)), ')');
                scatter(bdata(ind,1), bdata(ind,2), [], colors(ind,:), 'filled', markers(j), 'DisplayName', disp_name);
            end
        end
        text(bdata(:,1), bdata(:,2), olabel(:,1), 'VerticalAlignment','bottom','HorizontalAlignment','left');
    
    case 3,  % 2D plot, Multiplicity in size
        for i=1:length(unique_classes)
            ind = classes == unique_classes(i);
            scatter(bdata(ind,1), bdata(ind,2),sizes(ind), colors(ind,:),'filled','DisplayName',num2str(unique_classes(i)));
        end
        text(bdata(:,1), bdata(:,2), olabel(:,1), 'VerticalAlignment','bottom','HorizontalAlignment','left');
    
    case 4,  % 3D plot, Multiplicity in Z-axis
        for i=1:length(unique_classes)
            ind = classes == unique_classes(i);
            scatter3(bdata(ind,1), bdata(ind,2), mult(ind), [], colors(ind,:), 'filled', 'DisplayName',num2str(unique_classes(i)));
        end
        text(bdata(:,1), bdata(:,2), mult, olabel(:,1), 'VerticalAlignment','bottom','HorizontalAlignment','left');
    
    case 5,  % 3D plot, Multiplicity in size, classes in Z-axis
        for i=1:length(unique_classes)
            ind = classes == unique_classes(i);
            scatter3(bdata(ind,1), bdata(ind,2), normal_classes(ind), sizes(ind), colors(ind,:), 'filled', 'DisplayName',num2str(unique_classes(i)));
        end
        text(bdata(:,1), bdata(:,2), normal_classes, olabel(:,1), 'VerticalAlignment','bottom','HorizontalAlignment','left');
end

% Set axis labels and plot origin lines
xlabel(axlabel(1), 'FontSize', 16);
ylabel(axlabel(2), 'FontSize', 16);

ax = axis;
plot([0 0], ax(3:4), 'k--','HandleVisibility', 'off');
plot(ax(1:2), [0 0], 'k--','HandleVisibility', 'off');
axis(ax)
box on

  