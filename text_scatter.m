function text_scatter(fig_h,bdata,varargin)

% Print text in a Scatter plot.
%
% text_scatter(fig_h,bdata) % minimum call
% text_scatter(fig_h,bdata,'EleLabel',elabel,'ObsClass',classes,'Option',opt,'Multiplicity',mult,'BlurIndex',blur) % complete call
%
%
% INPUTS:
%
% fig_h: (1x1) figure handle
%
% bdata: (Nx2) bidimensional data 
%
% Optional INPUTS (parameter):
%
% 'EleLabel': [Nx1] name of the elements (numbers are used by default)
%
% 'ObsClass': [Nx1, str(N), {N}] groups for different visualization (a single
%   group by default)
%
% 'Option': (str or num) options for data plotting: binary code of the form 'ab' for:
%       a:
%           0: do not plot multiplicity
%           1: plot multiplicity
%       b: (for a 0)
%           0: filled marks
%           1: empty marks
%       b: (for a 1)
%           00: plot multiplicity info in the size of the markers.
%           01: plot multiplicity info in the form of the markers.
%           10: plot multiplicity information in the Z axis.
%           11: plot multiplicity info in the size of the markers and
%               classes in Z-axis
%
%   By deafult, opt = '00'. If less digits are specified, least significant
%   digits are set to 0, i.e. opt = 1 means a=1, b=00.
%
% 'Multiplicity': [Nx1] multiplicity of each row (1s by default)
%
% 'BlurIndex': [1x1] avoid blur when adding labels. The higher, the more labels
%   are printer (the higher blur). Inf shows all the labels (1 by default).
%
%
% OUTPUTS:
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
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
assert (nargin >= 2, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);

figure(fig_h);

N = size(bdata, 1);

% Introduce optional inputs as parameters (name-value pair) 
p = inputParser;
addParameter(p,'EleLabel',1:N);   
addParameter(p,'ObsClass',ones(N,1));   
addParameter(p,'Option','000');   
addParameter(p,'Multiplicity',ones(N,1)); 
addParameter(p,'BlurIndex',1);
parse(p,varargin{:});

% Extract inputs from inputParser for code legibility
elabel = p.Results.EleLabel;
opt = p.Results.Option;
classes = p.Results.ObsClass;
mult = p.Results.Multiplicity;
blur = p.Results.BlurIndex;

% Convert int arrays to str
if isnumeric(opt), opt=num2str(opt); end

% Convert row arrays to column arrays
if size(elabel,1)  == 1, elabel  = elabel';  end;
if size(classes,1) == 1, classes = classes'; end;
if size(mult,1) == 1, mult = mult'; end;

% Convert num arrays to str
if ~isempty(elabel) && isnumeric(elabel), elabel=num2str(elabel); end
if isnumeric(opt), opt=num2str(opt); end
if ~isempty(classes) && isnumeric(classes), classes=num2str(classes); end

% Complete opt
while length(opt)<3, opt = strcat(opt,'0'); end

% Convert char arrays to cell
if ischar(elabel),  elabel = cellstr(elabel); end;
if ischar(classes), classes = cellstr(classes); end;

% Validate dimensions of input data
assert(size(bdata,2) == 2, 'Dimension Error: paramter ''bdata'' must be N-by-2. Type ''help %s'' for more info.', routine(1).name);
if ~isempty(elabel), assert (isequal(size(elabel), [N 1]), 'Dimension Error: paramter ''EleLabel'' must be N-by-1. Type ''help %s'' for more info.', routine(1).name); end;
if ~isempty(classes), assert (isequal(size(classes), [N 1]), 'Dimension Error: parameter ''ObsClass'' must be N-by-1. Type ''help %s'' for more info.', routine(1).name); end;
assert (ischar(opt) && length(opt)==3, 'Dimension Error: parameter ''Option'' must be a string or num of maximum 3 bits. Type ''help %s'' for more info.', routine(1).name);
if ~isempty(mult), assert (isequal(size(mult), [N 1]), 'Dimension Error: parameter ''Multiplicity'' must be N-by-1. Type ''help %s'' for more info.', routine(1).name); end;
if ~isempty(blur), assert (isequal(size(blur), [1 1]), 'Dimension Error: parameter ''BlurIndex'' must be 1-by-1. Type ''help %s'' for more info.', routine(1).name); end;

% Validate values of input data
assert (isempty(find(opt~='0' & opt~='1')), 'Value Error: parameter ''Option'' must contain binary values. Type ''help %s'' for more info.', routine(1).name);


%% Main code

% Get ordering of classes
unique_classes = unique(classes,'stable');
if iscell(classes)
     ord_classes = arrayfun(@(x) find(strcmp(unique_classes, x), 1), classes);
else
     ord_classes = arrayfun(@(x) find(unique_classes == x, 1), classes);
end
unique_ord_classes = unique(ord_classes);

ax = axis;
deltax = (ax(2)-ax(1))/100;
deltay = (ax(4)-ax(3))/100;

bdata = bdata + 0.01*randn(size(bdata)).*(ones(size(bdata,1),1)*std(bdata));
mar = 0.1;
if ~isempty(elabel)
    for i=1:N
        suffx = length(char(strtrim(elabel(i,1))));
        ind = [1:(i-1) (i+1):size(bdata,1)];

        dx = (bdata(ind,1)-bdata(i,1))/deltax;
        dxM = dx;
        dxM(dxM<0) = Inf;
        dxm = dx;
        dxm(dxm>mar) = Inf;
        dxm(dxm>0 & dxm<=mar) = 0;
        dy = (bdata(ind,2)-bdata(i,2))/deltay;
        dyM = dy;
        dyM(dyM<0) = Inf;
        dym = dy;
        dym(dym>mar) = Inf;
        dym(dym>0 & dym<=mar) = 0;

        % Labels in any direction: not used

%         d = min([dxM.^2+dyM.^2 dxM.^2+dym.^2 dxm.^2+dyM.^2 dxm.^2+dym.^2]);
%         if length(find(d > 10/blur))>1 || isempty(ind),
%             quad = find(d==max(d),1);
%             switch quad,
%                 case 1,
%                     posx = bdata(i,1)+deltax;
%                     posy = bdata(i,2)+deltay;
%                 case 2,
%                     posx = bdata(i,1)+deltax;
%                     posy = bdata(i,2)-6*deltay;
%                 case 3,
%                     posx = bdata(i,1)-deltax-suffx/2;
%                     posy = bdata(i,2)+2*deltay;
%                 case 4,
%                     posx = bdata(i,1)-deltax-suffx/2;
%                     posy = bdata(i,2)-6*deltay;
%             end
%
%
%         % Labels only to the right: used
%
        d = min([dxM.^2+dyM.^2 dxM.^2+dym.^2 dxm.^2+dyM.^2 dxm.^2+dym.^2]);
        if (length(find(d > 10/blur))>1 && length(find(d(1:2) > 10/blur))>0)|| isempty(ind)
            quad = find(d(1:2)==max(d(1:2)),1);
            switch quad
                case 1
                    posx = bdata(i,1)+deltax;
                    posy = bdata(i,2)+deltay;
                case 2
                    posx = bdata(i,1)+deltax;
                    posy = bdata(i,2)-6*deltay;
                case 3
                    posx = bdata(i,1)-deltax-suffx/2;
                    posy = bdata(i,2)+2*deltay;
                case 4
                    posx = bdata(i,1)-deltax-suffx/2;
                    posy = bdata(i,2)-6*deltay;
            end

%         % Labels only to upper right: not used
%
%         d = min([dxM.^2+dyM.^2 dxM.^2+dym.^2 dxm.^2+dyM.^2 dxm.^2+dym.^2]);
%         if (length(find(d > 10/blur))>1 && d(1) > 10/blur)|| isempty(ind),
%             posx = bdata(i,1)+deltax;
%             posy = bdata(i,2)+deltay;


            switch opt
                case '110'
                    text(posx, posy, mult(i), strtrim(elabel(i,1)),'VerticalAlignment','bottom', 'HorizontalAlignment','left','FontSize', 12);
                case '111'
                    text(posx, posy, ord_classes(i), strtrim(elabel(i,1)),'VerticalAlignment','bottom', 'HorizontalAlignment','left','FontSize', 12);
                otherwise
                    text(posx, posy, strtrim(elabel(i,1)),'VerticalAlignment','bottom', 'HorizontalAlignment','left','FontSize', 12);
            end
        end
    end
end
