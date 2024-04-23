function fig_h = plot_vec(vec,varargin)

% Bar or line plot.
%
% plot_vec(vec) % minimum call
% plot_vec(vec,'PARAM1',val1,'PARAM2',val2,...)
%
%
% INPUTS:
%
% vec: [NxM] vector/s to plot. 
%
% Optional INPUTS (Parameters):
%
% 'EleLabel': [Nx1] name of the vector elements (numbers are used by default)
%
% 'ObsClass': [Nx1, str(N), {N}] groups for different visualization (a single 
%   group by default)
%
% 'XYLabel': {2} xlabel and ylabel (nothing by default)
%
% 'LimCont': [NxL or Lx1] L control limits (nothing by default)
%
% 'Option': (str or num) options for data plotting: binary code of the form 'ab' for:
%       a:
%           0: line plot
%           1: bar plot
%       b: 
%           0: plot for numerical classes (consistent with a colorbar)
%           1: plot for categorical classes (consistent with a legend)
%
%   By deafult, opt = '11'. If less digits are specified, least significant
%   digits are set to 0, i.e. opt = 1 means a=1, b=0
%
% 'VecLabel': [Mx1] name of the vectors (numbers are used by default)
%
% 'Multiplicity': [NxM] multiplicity of each row (1s by default)
%
% 'Markers': [1x3] thresholds for the different marker size (20, 50 and 100 by default)
%
%
% OUTPUTS:
%
% fig_h: (1x1) figure handle.
%
%
% EXAMPLE OF USE: To plot three lines with constant control limits:
%
% fig_h = plot_vec(randn(100,3),'XYLabel',{'Functions','Time'},'LimCont',[1, -1, 3]);
%
%
% EXAMPLE OF USE: with labels and classes in observations and variable limit:
%
% fig_h = plot_vec(randn(5,3),'EleLabel',{'one','two','three','four','five'},'ObsClass',[1 1 1 2 2],'XYLabel',{[],'Functions'},'LimCont',randn(5,1),'Option','11');
%
%
% EXAMPLE OF USE: with labels, multiplicity and classes in observations and variable limit:
%
% fig_h = plot_vec(randn(5,3),'EleLabel',{'one','two','three','four','five'},'ObsClass',[1 1 1 2 2],'XYLabel',{[],'Functions'},'LimCont',randn(5,1),'Option',11,'Multiplicity',100*rand(5,1),'Markers',[20 50 100]);
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
%           Alejandro Perez Villegas (alextoni@gmail.com)
% last modification: 23/Apr/2024
%
% Copyright (C) 2024  University of Granada, Granada
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
if size(vec,1) == 1,     vec = vec'; end;
N = size(vec, 1);
M = size(vec, 2);

% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
addParameter(p,'EleLabel',1:N);   
addParameter(p,'ObsClass',[]);
addParameter(p,'XYLabel',{'',''});
addParameter(p,'LimCont',[]);
addParameter(p,'Option','11');
addParameter(p,'Multiplicity',ones(N,1));
addParameter(p,'Markers',[20 50 100]);
addParameter(p,'VecLabel',1:M);
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
elabel = p.Results.EleLabel;
classes = p.Results.ObsClass;
xylabel = p.Results.XYLabel;
lcont = p.Results.LimCont;
opt = p.Results.Option;
mult = p.Results.Multiplicity;
maxv = p.Results.Markers;
vlabel = p.Results.VecLabel;

% Convert num arrays to str
if isnumeric(opt), opt=num2str(opt); end

% Correct for opt integrity
while length(opt)<2, opt = strcat(opt,'0'); end
if opt(2) == 0 && ~isnumeric(classes), opt(2) = 1; end

% Convert row arrays to column arrays
if size(elabel,1)  == 1, elabel = elabel'; end;
if size(classes,1) == 1, classes = classes'; end;
if size(lcont,1) == 1, lcont = lcont'; end;
if size(vlabel,1)  == 1, vlabel = vlabel'; end;
if size(mult,1) == 1, mult = mult'; end;
if size(maxv,2) == 1, maxv = maxv'; end;

% Convert num arrays to str
if ~isempty(vlabel) && isnumeric(vlabel), vlabel=num2str(vlabel); end
if ~isempty(classes) && isnumeric(classes) && opt(2)=='1', classes=num2str(classes); end

% Convert char arrays to cell
if ischar(elabel),  elabel = cellstr(elabel); end;
if ischar(classes), classes = cellstr(classes); end;
if ischar(xylabel),  xylabel = cellstr(xylabel); end;
if ischar(vlabel),  vlabel = cellstr(vlabel); end;

% Validate dimensions of input data

if ~isempty(elabel), assert (isequal(size(elabel), [N 1]), 'Dimension Error: parameter ''EleLabel'' must be N-by-1. Type ''help %s'' for more info.', routine(1).name); end;
if ~isempty(classes), assert (isequal(size(classes), [N 1]), 'Dimension Error: parameter ''ObsClass'' must be N-by-1. Type ''help %s'' for more info.', routine(1).name); end;
if ~isempty(xylabel), assert (length(xylabel) == 2, 'Dimension Error: parameter ''XYLabel'' must contain 2 cell elements. Type ''help %s'' for more info.', routine(1).name); end;
if ~isempty(lcont), assert (isequal(size(lcont,1), N) || isequal(size(lcont,2), 1), 'Dimension Error: parameter ''LimCont'' must be N-by-L or L-by-1. Type ''help %s'' for more info.', routine(1).name); end;
assert (ischar(opt) && length(opt)==2, 'Dimension Error: parameter ''Option'' must be a string or num of maximum 2 bits. Type ''help %s'' for more info.', routine(1).name);
if ~isempty(vlabel), assert (isequal(size(vlabel), [M 1]), 'Dimension Error: parameter ''VecLabel'' must be M-by-1. Type ''help %s'' for more info.', routine(1).name); end;
if ~isempty(mult), assert (isequal(size(mult), [N 1]), 'Dimension Error: parameter ''Multiplicity'' must be N-by-1. Type ''help %s'' for more info.', routine(1).name); end;
if ~isempty(maxv), assert (isequal(size(maxv), [1 3]), 'Dimension Error: parameter ''Markers'' must be 1-by-3. Type ''help %s'' for more info.', routine(1).name); end;
 
% Validate values of input data
assert (isempty(find(opt~='0' & opt~='1')), 'Value Error: parameter ''Option'' must contain binary values. Type ''help %s'' for more info.', routine(1).name);
   
% Convert constant limits in vectors
if ~isempty(lcont) && ~isequal(size(lcont,1), N), lcont = (lcont*ones(1,N))'; end;
    
% Exception: bar plot with multivariate vec and one-observation class  
if ~opt(1) && ~isempty(classes) && size(vec, 2)>1
    unique_classes = unique(classes);
    assert (min(hist(classes,unique(classes)))>1, 'Exception: Cannot visualize a multivariate bar plot with one-observation classes. Try setting the 6th argument to 1.'); 
end;

%% Main code

% Create figure window
fig_h = figure;
hold on;

% Sort data for colorbar
if opt(2)=='0'
    [classes,ord] = sort(classes,'ascend');
    cax = [min(classes) max(classes)];
    classes = num2str(classes); 
    classes = cellstr(classes);
    vec = vec(ord,:);
    elabel = elabel(ord);
    mult = mult(ord);
end

% Get ordering of classes
unique_classes = unique(classes,'stable');
if iscell(classes)
    ord_classes = arrayfun(@(x) find(strcmp(unique_classes, x), 1), classes);
else
    ord_classes = arrayfun(@(x) find(unique_classes == x, 1), classes);
end
unique_ord_classes = unique(ord_classes);

% Plot vectors

bins = [0 1 maxv Inf];

sizes = [];
for i=1:length(bins)-1
    sizes (i) = round(.5 * i^2 * pi);
end

% Plot multiplicity
for j=1:length(bins)-1
    ind = mult>bins(j) & mult<=bins(j+1);
    if isnumeric(elabel)
        plot(elabel(find(ind)), 0*find(ind), 'ko', 'MarkerSize', sizes(j), 'HandleVisibility', 'off');
    else
        plot(find(ind), 0*find(ind), 'ko', 'MarkerSize', sizes(j), 'HandleVisibility', 'off');
    end
end

if ~isempty(classes)
    if opt(2) == '0'
        color_list = parula(length(unique_ord_classes));
    else
        color_list = hsv(length(unique_ord_classes));
    end
    for i=1:length(unique_ord_classes)
        ind = ord_classes == unique_ord_classes(i);
        if isnumeric(elabel) && length(elabel)==length(unique(elabel))
            vind = elabel(find(ind));
        else
            vind = find(ind);
        end
            
        if opt(1) == '0'
           plot(vind, vec(ind,:), 'Color', 'none', 'Marker','O', 'MarkerFaceColor', color_list(i,:), 'DisplayName', unique_classes{i});
        else 
           bar(vind, vec(ind,:), 0.8, 'FaceColor', color_list(i,:), 'EdgeColor', 'none', 'DisplayName', unique_classes{i});
        end
    end 
else
    color_list = hsv(M);
    for i=1:M
        if opt(1) == '0'
            if isnumeric(elabel) && length(elabel)==length(unique(elabel))
                plot(elabel, vec(:,i), 'LineWidth', 2, 'Color', color_list(i,:), 'DisplayName', vlabel{i});
            else
                plot(vec(:,i), 'LineWidth', 2, 'Color', color_list(i,:), 'DisplayName', vlabel{i});
            end
        else
            if isnumeric(elabel) && length(elabel)==length(unique(elabel))
                bar(elabel, vec(:,i), 'FaceColor', color_list(i,:), 'EdgeColor', 'none', 'DisplayName', vlabel{i});
            else
                bar(vec(:,i), 'FaceColor', color_list(i,:), 'EdgeColor', 'none', 'DisplayName', vlabel{i});
            end
        end
    end    
end

% Plot control limits
if ~isempty(lcont)
    hold on
    b = [0.5:(N+1);0.5:(N+1)];
    for i=1:size(lcont,2)
        a = [lcont(:,i)';lcont(:,i)'];
        plot(b(2:(end-1))',a(:),'r--','LineWidth',2,'HandleVisibility', 'off');
    end
end    

% Get axes handler
axes_h = get(fig_h,'Children');
for i=1:length(axes_h)
    if strcmp(get(axes_h(i), 'type'), 'axes')
        set(axes_h(i), 'FontSize', 14);
        val=i;
    end
end
axes_h = axes_h(i); 

% Set ticks and labels
if ~isempty(elabel) & ~isnumeric(elabel),
    label_length = max(cellfun('length', elabel));
    label_size = 300/(length(find(~cellfun('isempty', elabel)))*label_length);
    set(axes_h, 'FontSize', max(min(14,round(label_size)), 10));
    stepN = ceil(0.2*N/label_size);
    if stepN==1,
        vals = 1:N;
        set(axes_h,'XTick',vals);
        set(axes_h,'XTickLabel',elabel(vals));
    else
        set(axes_h,'XTickMode','auto');
        set(axes_h, 'FontSize', 14);
    end
end

if ~isempty(xylabel)
    xlabel(xylabel{1}, 'FontSize', 16);
    ylabel(xylabel{2}, 'FontSize', 16);
end

% Set axis
axis tight;
ax = axis;
axis auto;
ax2 = axis;
axis([ax(1:2) ax2(3:4)]);

%legend off
box on;
hold off;
        