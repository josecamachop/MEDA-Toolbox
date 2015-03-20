
function fig_h = plot_vec2 (vec,olabel,classes,slabel,lcont)

% Bar plot.
%
% plot_vec(vec) % minimum call
% plot_vec(vec,olabel,classes,slabel,lcont,fig_h) % complete call
%
%
% INPUTS:
%
% vec: (Nx1) vector to plot. 
%
% olabel: (Nx1) name of the x-variables
%
% classes: (Nx1) grouping variable used to classify observations/variables
%
% slabel: (str) y-variable/statistic plotted (nothing by default)
%
% lcont: (2xN) control limits.
%
%
% OUTPUTS:
%
% fig_h: (1x1) figure handle.
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
%           Alejandro Perez Villegas (alextoni@gmail.com)
% last modification: 19/Mar/2015.
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
assert (nargin >= 1, 'Error: Missing arguments.');
if nargin < 2 || isempty(olabel), olabel = []; end;
if nargin < 3 || isempty(classes), classes = []; end;
if nargin < 4 || isempty(slabel), slabel = ''; end;
if nargin < 5 || isempty(lcont),  lcont = []; end;

% Convert char arrays to cell
if ischar(olabel),  olabel = cellstr(olabel); end;
if ischar(classes), classes = cellstr(classes); end;
if ischar(slabel),  slabel = cellstr(slabel); end;

% Convert row arrays to column arrays
if size(vec,1) == 1,     vec = vec'; end;
if size(olabel,1)  == 1, olabel = olabel'; end;
if size(classes,1) == 1, classes = classes'; end;

% Validate dimensions of input data
N = size(vec, 1);
assert (isequal(size(vec), [N 1]), 'Dimension Error: vec must be N-by-1.');
if ~isempty(olabel),  assert (isequal(size(olabel), [N 1]),  'Dimension Error: label must be N-by-1.'); end;
if ~isempty(classes), assert (isequal(size(classes), [N 1]), 'Dimension Error: classes must be N-by-1.'); end;
if ~isempty(slabel),  assert (isequal(size(slabel), [1 1]),  'Dimension Error: slabel must be a single string.'); end;
if ~isempty(lcont),   assert (isequal(size(lcont), [2 N]),   'Dimension Error: lcont must be 2-by-N.'); end;



%% Main code

% Create figure window
fig_h = figure;
hold on;


% Plot bars
if ~isempty(classes)
    unique_classes = unique(classes);
    color_list = hsv(length(unique_classes));
    for i=1:length(unique_classes)
        ind = classes == unique_classes(i);
        bar(find(ind), vec(ind,:), 'FaceColor', color_list(i,:), 'EdgeColor', 'none', 'DisplayName', num2str(unique_classes(i)));
    end
else
    bar(vec, 'r', 'EdgeColor', 'none');
end


% Plot control limits
if ~isempty(lcont)
    hold on
    b = [0.5:(N+1);0.5:(N+1)];
    a = [lcont(1,:);lcont(1,:)];
    plot(b(2:(end-1))',a(:),'r--','LineWidth',2,'HandleVisibility', 'off');
    a = [lcont(2,:);lcont(2,:)];
    plot(b(2:(end-1))',a(:),'r--','LineWidth',2,'HandleVisibility', 'off');
end    

% Get axes handler
axes_h = get(fig_h,'Children');
if length(axes_h)>1, axes_h = axes_h(1); end;

% Set ticks and labels
label_size = max(min(14,round(300/length(olabel))), 9);
set(axes_h, 'FontSize', label_size);
if ~isempty(olabel)
    set(axes_h,'XTick',1:N);
    label_length = max(cellfun('length', olabel));
    if label_length > 2
        % rotate labels
        set(axes_h,'XTickLabel',[]);
        b = get(axes_h, 'XTick');
        c = get(axes_h, 'YTick');
        text(b, repmat(c(1)-.1*(c(2)-c(1)),length(b),1), olabel,...
            'FontSize',label_size,'HorizontalAlignment','right','rotation',90);
    else
        % dont rotate labels
        set(axes_h,'XTickLabel', olabel);
    end
end
if ~isempty(slabel)
    ylabel(slabel, 'FontSize', 16);
end

% Set axis
axis tight
ax = axis;
axis auto
ax2 = axis;
axis([ax(1:2) ax2(3:4)])

       