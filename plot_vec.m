
function fig_h = plot_vec(vec,elabel,classes,yxlabel,lcont,opt,vlabel)

% Bar plot.
%
% fig_h = plot_vec(vec) % minimum call
% fig_h = plot_vec(vec,elabel,classes,yxlabel,lcont,opt,vlabel) % complete call
%
%
% INPUTS:
%
% vec: [NxM] vector/s to plot. 
%
% elabel: [Nx1] name of the vector elements (numbers are used by default)
%
% classes: [Nx1, str(N), {N}] groups for different visualization (a single 
%   group by default)
%
% yxlabel: {1 or 2} ylabel and xlabel. If only one label is specified, it 
%   is understood as the ylabel (nothing by default)
%
% lcont: [NxL or Lx1] L control limits (nothing by default)
%
% opt: [1x1] options for data plotting
%       0: bar plot (by default)
%       otherwise: line plot
%
% vlabel: [Mx1] name of the vectors (numbers are used by default)
%
%
% OUTPUTS:
%
% fig_h: [1x1] figure handle
%
%
% EXAMPLE OF USE: To plot three lines with constant control limits:
%
% fig_h = plot_vec(randn(100,3),[],[],{'Functions','Time'},[1, -1, 3], 1);
%
%
% EXAMPLE OF USE: with labels and classes in observations and variable limit:
%
% fig_h = plot_vec(randn(5,3),{'one','two','three','four','five'},[1 1 1 2 2],{'Functions'},randn(5,1));
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
%           Alejandro Perez Villegas (alextoni@gmail.com)
% last modification: 30/Mar/2016.
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
assert (nargin >= 1, 'Error in the number of arguments. Type ''help %s'' for more info.', routine.name);
if size(vec,1) == 1,     vec = vec'; end;
N = size(vec, 1);
M = size(vec, 2);
if nargin < 2 || isempty(elabel), elabel = 1:N; end;
if nargin < 3 || isempty(classes), classes = []; end;
if nargin < 4 || isempty(yxlabel), yxlabel = ''; end;
if nargin < 5 || isempty(lcont),  lcont = []; end;
if nargin < 6 || isempty(opt),  opt = 0; end;
if nargin < 7 || isempty(vlabel),  vlabel = 1:M; end;

% Convert row arrays to column arrays
if size(elabel,1)  == 1, elabel = elabel'; end;
if size(classes,1) == 1, classes = classes'; end;
if size(lcont,1) == 1, lcont = lcont'; end;
if size(vlabel,1)  == 1, vlabel = vlabel'; end;

% Convert int arrays to str
if ~isempty(elabel) && isnumeric(elabel), elabel=num2str(elabel); end
if ~isempty(vlabel) && isnumeric(vlabel), vlabel=num2str(vlabel); end

% Convert char arrays to cell
if ischar(elabel),  elabel = cellstr(elabel); end;
if ischar(classes), classes = cellstr(classes); end;
if ischar(yxlabel),  yxlabel = cellstr(yxlabel); end;
if ischar(vlabel),  vlabel = cellstr(vlabel); end;

% Validate dimensions of input data
if ~isempty(elabel), assert (isequal(size(elabel), [N 1]), 'Dimension Error: 2nd argument must be N-by-1. Type ''help %s'' for more info.', routine.name); end;
if ~isempty(classes), assert (isequal(size(classes), [N 1]), 'Dimension Error: 3rd argument must be N-by-1. Type ''help %s'' for more info.', routine.name); end;
if ~isempty(yxlabel), assert (length(yxlabel) <= 2, 'Dimension Error: 4th argument must contain 2 cell elements at most. Type ''help %s'' for more info.', routine.name); end;
if ~isempty(lcont), assert (isequal(size(lcont,1), N) | isequal(size(lcont,2), 1), 'Dimension Error: 5th argument must be N-by-L or L-by-1. Type ''help %s'' for more info.', routine.name); end;
assert (isequal(size(opt), [1 1]), 'Dimension Error: 6th argument must be 1-by-1. Type ''help %s'' for more info.', routine.name);
if ~isempty(vlabel), assert (isequal(size(vlabel), [M 1]), 'Dimension Error: 7th argument must be M-by-1. Type ''help %s'' for more info.', routine.name); end;
    
% Convert constant limits in vectors
if ~isempty(lcont) && ~isequal(size(lcont,1), N), lcont = (lcont*ones(1,N))'; end;
    
% Exception: bar plot with multivariate vec and one-observation class  
if ~opt && ~isempty(classes) && size(vec, 2)>1,
    unique_classes = unique(classes);
    assert (min(hist(classes,unique(classes)))>1, 'Exception: Cannot visualize a multivariate bar plot with one-observation classes. Try setting the 6th argument to 1.'); 
end;
    

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

% Plot vectors

if ~isempty(classes)
    unique_classes = unique(classes);
    color_list = hsv(length(unique_classes));
    if opt,
        plot(vec,'k','HandleVisibility', 'off');
    end  
    for i=1:length(unique_classes)
        ind = classes == unique_classes(i);
        if opt,
            plot(find(ind), vec(ind,:), 'Color', 'none', 'Marker','O', 'MarkerFaceColor', color_list(i,:), 'DisplayName', num2str(unique_classes(i)));
        else 
            bar(find(ind), vec(ind,:), 'FaceColor', color_list(i,:), 'EdgeColor', 'none', 'DisplayName', num2str(unique_classes(i)));
        end
    end
else
    color_list = hsv(M);
    for i=1:M,
        if opt,
            plot(vec(:,i), 'Color', color_list(i,:), 'DisplayName', vlabel{i});
        else
            bar(vec(:,i), 'FaceColor', color_list(i,:), 'EdgeColor', 'none', 'DisplayName', vlabel{i});
        end
    end
end

% Plot control limits
if ~isempty(lcont)
    hold on
    b = [0.5:(N+1);0.5:(N+1)];
    for i=1:size(lcont,2),
        a = [lcont(:,i)';lcont(:,i)'];
        plot(b(2:(end-1))',a(:),'r--','LineWidth',2,'HandleVisibility', 'off');
    end
end    

% Get axes handler
axes_h = get(fig_h,'Children');
if length(axes_h)>1, axes_h = axes_h(1); end;

% Set ticks and labels
label_size = max(min(14,round(300/length(elabel))), 9);
set(axes_h, 'FontSize', label_size);
if ~isempty(elabel)
    set(axes_h,'XTick',1:N);
    label_length = max(cellfun('length', elabel));
    
    set(axes_h,'XTickLabel',elabel);
end
if ~isempty(yxlabel)
    ylabel(yxlabel{1}, 'FontSize', 16);
    if length(yxlabel)>1,
        xlabel(yxlabel{2}, 'FontSize', 16);
    end
end

% Set axis
axis tight
ax = axis;
axis auto
ax2 = axis;
axis([ax(1:2) ax2(3:4)])

legend off
box on
hold off

       