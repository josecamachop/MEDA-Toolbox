
function fig_h = plot_scatter(bdata,elabel,classes,xylabel,lcont,opt)

% Scatter plot.
%
% fig_h = plot_scatter(bdata) % minimum call
% fig_h = plot_scatter(bdata,elabel,classes,xylabel,lcont,opt) % complete call
%
%
% INPUTS:
%
% bdata: [Nx2] bidimensional data to plot. 
%
% elabel: [Nx1] name of the elements (numbers are used by default)
%
% classes: [Nx1, str(N), {N}] groups for different visualization (a single 
%   group by default)
%
% xylabel: {2} xlabel and ylabel (nothing by default)
%
% lcont: {2} control limits on x and y axis (nothing by default)
%
% opt: [1x1] options for data plotting
%       0: filled marks (by default)
%       1: empty marks
%
%
% OUTPUTS:
%
% fig_h: (1x1) figure handle
%
%
% EXAMPLE OF USE: Plot random data with empty marks and control limits:
%
% fig_h = plot_scatter(rand(100,2),[],[],{'Y','X'},{0.8,0.8},1);
%
%
% EXAMPLE OF USE: with labels and classes in elements:
%
% fig_h = plot_scatter(randn(5,2),{'one','two','three','four','five'},[1 1 1 2 2],{'Y','X'});
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
%           Alejandro Perez Villegas (alextoni@gmail.com)
% last modification: 31/Mar/2016.
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

% Set default values
routine=dbstack;
assert (nargin >= 1, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
N = size(bdata, 1);
if nargin < 2 || isempty(elabel), elabel = 1:N; end;
if nargin < 3 || isempty(classes), classes = []; end;
if nargin < 4 || isempty(xylabel), xylabel = ''; end;
if nargin < 5 || isempty(lcont),  lcont = []; end;
if nargin < 6 || isempty(opt),  opt = 0; end;

% Convert row arrays to column arrays
if size(elabel,1)  == 1, elabel = elabel'; end;
if size(classes,1) == 1, classes = classes'; end;
if size(lcont,1) == 1, lcont = lcont'; end;

% Convert int arrays to str
if ~isempty(elabel) && isnumeric(elabel), elabel=num2str(elabel); end

% Convert char arrays to cell
if ischar(elabel),  elabel = cellstr(elabel); end;
if ischar(classes), classes = cellstr(classes); end;
if ischar(xylabel),  xylabel = cellstr(xylabel); end;

% Validate dimensions of input data
if ~isempty(elabel), assert (isequal(size(elabel), [N 1]), 'Dimension Error: 2nd argument must be N-by-1. Type ''help %s'' for more info.', routine(1).name); end;
if ~isempty(classes), assert (isequal(size(classes), [N 1]), 'Dimension Error: 3rd argument must be N-by-1. Type ''help %s'' for more info.', routine(1).name); end;
if ~isempty(xylabel), assert (length(xylabel) == 2, 'Dimension Error: 4th argument must contain 2 cell elements. Type ''help %s'' for more info.', routine(1).name); end;
if ~isempty(lcont), assert (iscell(lcont) & isequal(size(lcont), [2 1]), 'Dimension Error: 5th argument must be a cell of 2 elements. Type ''help %s'' for more info.', routine(1).name); end;
assert (isequal(size(opt), [1 1]), 'Dimension Error: 6th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
    

%% Main code

% Create figure window
fig_h = figure;
hold on;

% Preprocess classes to force them start with 1, 2...n,
unique_classes = unique(classes);
if iscell(classes)
    classes = arrayfun(@(x) find(strcmp(unique_classes, x), 1), classes);
else
    classes = arrayfun(@(x) find(unique_classes == x, 1), classes);
end

% Plot points
a = gscatter(bdata(:,1), bdata(:,2), classes, [], 'o');

% Fill marks
if opt == 0
    for i=1:length(a)
        color = get(a(i), 'Color');
        set(a(i), 'MarkerFaceColor',color);
    end
end

% Plot labels
ax = axis;
f = 5;
deltax = (ax(2)-ax(1))/100;
deltay = (ax(4)-ax(3))/100;
if ~isempty(elabel)
    for i=1:N
        nch = length(char(strtrim(elabel(i,1))));
        if length(find((bdata(:,1)>bdata(i,1)) & (bdata(:,1)<bdata(i,1)+deltax*nch*f) & (bdata(:,2)<bdata(i,2)+deltay*f) & (bdata(:,2)>bdata(i,2)-deltay*f)))<1,
            text(bdata(i,1)+deltax, bdata(i,2)+deltay, strtrim(elabel(i,1)),'VerticalAlignment','bottom', 'HorizontalAlignment','left');
        end
    end
end

ax = axis;
ax([1 3]) = min(ax([1 3]),zeros(1,2));
ax([2 4]) = max(ax([2 4]),zeros(1,2));

if ~isempty(lcont) % Plot control limits
    if ~isempty(lcont{2})
        ax(3) = min([ax(3);lcont{2}(:)]);
        ax(4) = max([ax(4);lcont{2}(:)]);
    end
    if ~isempty(lcont{1})
        ax(1) = min([ax(1);lcont{1}(:)]);
        ax(2) = max([ax(2);lcont{1}(:)]);
        for i=1:length(lcont{1}),
            plot([lcont{1}(i) lcont{1}(i)], ax(3:4), 'r--','LineWidth',2, 'HandleVisibility', 'off');
        end
    end
    if ~isempty(lcont{2})
        for i=1:length(lcont{2}),
            plot(ax(1:2),[lcont{2}(i) lcont{2}(i)], 'r--','LineWidth',2, 'HandleVisibility', 'off');
        end
    end    
else % Plot origin lines
    plot([0 0], ax(3:4), 'k--', 'HandleVisibility', 'off');
    plot(ax(1:2), [0 0], 'k--', 'HandleVisibility', 'off');
end    
axis(ax)

% Set axis labels
if ~isempty(xylabel)
    xlabel(xylabel{1}, 'FontSize', 16);
    ylabel(xylabel{2}, 'FontSize', 16);
end



legend off
box on
hold off

