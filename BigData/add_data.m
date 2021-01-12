function add_data(name,path,data,label,class,type,thres,preci,debug)

% Add data to a file in the clustering file system. 
%
% add_data(name,path,data,label,class,type,thres) % minimum call
% add_data(name,path,data,label,class,type,thres,preci,debug) % complete call
%
%
% INPUTS:
%
% name: (str) name of the file.
%
% path: (str) path to the directory where the clustering data files are
%   located.
%
% data: [NxM] observations to include in the file.
%
% label: [Nx1] name of the observations (filenames are used by default)
%
% class: [1x1] class associated to the observations.
%
% type: [1x1] type of update:
%       'w'     open file for writing; discard existing contents
%       'a'     open or create file for writing; append data to end of file
%
% thres: [1x1] maximum number of entries in a file.
%
% preci: [1x1] number of decimals (8 by default)
%
% debug: [1x1] disply debug messages
%       0: no messages are displayed.
%       1: display only main messages (default). In the present routine, no 
%           messages are displayed.
%       2: display all messages.
%
%
% OUTPUTS:
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 12/Jan/2021
%
% Copyright (C) 2021  University of Granada, Granada
% Copyright (C) 2021  Jose Camacho Paez
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
    
%% Arguments checking

% Set default values
routine=dbstack;
assert (nargin >= 7, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
N = size(data, 1);
M = size(data, 2);
if nargin < 8 || isempty(preci), preci=8; end;
if nargin < 9 || isempty(debug), debug = 1; end;

% Convert row arrays to column arrays
if size(label,1)  == 1, label = label'; end;

% Validate dimensions of input data
assert (isequal(size(label), [N 1]), 'Dimension Error: 4th argument must be N-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(class), [1 1]), 'Dimension Error: 5th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(type), [1 1]), 'Dimension Error: 6th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(thres), [1 1]), 'Dimension Error: 7th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(preci), [1 1]), 'Dimension Error: 8th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(size(debug), [1 1]), 'Dimension Error: 9th argument must be 1-by-1. Type ''help %s'' for more info.', routine(1).name);

% Validate values of input data
assert (type=='w' || type=='a', 'Value Error: 6th argument must be ''w'' or ''a''. Type ''help %s'' for more info.', routine(1).name);
assert (thres>0, 'Value Error: 7th argument must be above 0. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(fix(thres), thres), 'Value Error: 7th argument must contain an integer. Type ''help %s'' for more info.', routine(1).name);
assert (preci>0, 'Value Error: 8th argument must be above 0. Type ''help %s'' for more info.', routine(1).name);
assert (isequal(fix(preci), preci), 'Value Error: 8th argument must contain an integer. Type ''help %s'' for more info.', routine(1).name);
assert (debug==0 || debug==1 || debug==2, 'Value Error: 9th argument must be 0, 1 or 2. Type ''help %s'' for more info.', routine(1).name);


%% Main code

preci_str = sprintf('%%.%df,',preci);

s=size(data);
file=[path name '.txt'];

if debug>1, disp(['add data in file: ' file ' ...']), end;

if isequal('a',type),
    fid=fopen(file,'r');
    a = fscanf(fid,'%d',3);
    lev=a(1);
    s2=a(2);
    if lev == 0,
        fclose(fid);
        stot = s(1) + s2;
        if stot > thres,
            [data2, label2] = read_data(name,path,s(2),debug);
            data = [data2;data];
            label = {label{:} label2{:}};
            add_data1(name,path,data,label,class,'w',thres,1,preci);
        else
            fid=fopen(file,'r+');
            str=sprintf('%d %d %d',0,stot,class); 
            str = [str char(12*ones(1,10-length(str)))];
            fprintf(fid,'%s\n',str);  
            fseek(fid,0,'eof');
            for u=1:s(1),
                a=num2str(data(u,:),preci_str);
                i=find(~isspace(a));
                a=a(i);
                fprintf(fid,'%s: %s\n',label{u},a);
            end
            fclose(fid);
        end
    else
        for i=1:s2,
            name2 = fscanf(fid,'%s',1);
        end  
        fclose(fid);
        [data2, label2] = read_data(name2,path,s(2),debug);
        data = [data2;data];
        label = {label{:} label2{:}};
        add_data1(name,path,data,label,class,'a',thres,s2,preci);
    end
else
    if s(1) > thres,
        add_data1(name,path,data,label,class,'w',thres,1,preci);
    else
        fid=fopen(file,'w');
        str=sprintf('%d %d %d',0,s(1),class); 
        str = [str char(12*ones(1,10-length(str)))];
        fprintf(fid,'%s\n',str);  
        for u=1:s(1),
            a=num2str(data(u,:),preci_str);
            i=find(~isspace(a));
            a=a(i);
            fprintf(fid,'%s: %s\n',label{u},a);
        end
        fclose(fid);
    end
end


function add_data1(name,path,data,label,class,type,thres,s2,preci)

if nargin < 8, error('Error in the number of arguments.'); end;
if nargin < 9, preci=8; end;

preci_str = sprintf('%%.%df,',preci);

s=size(data);   
nfich = ceil(s(1)/thres);
file=[path name '.txt'];
if isequal('a',type),
    fid=fopen(file,'r+');
else
    fid=fopen(file,'w');
end
str=sprintf('%d %d %d',1,nfich+s2-1,class);
str = [str char(12*ones(1,10-length(str)))];
fprintf(fid,'%s\n',str);
fseek(fid,0,'eof');


for i=s2:nfich+s2-1,
    if ~(i==s2 && isequal('a',type)), fprintf(fid,'%s_%d\n',name,i); end;
    file=[path name '_' num2str(i) '.txt'];
    fid2=fopen(file,'w');
    indu = (((i-s2)*thres+1):min(s,(i-s2+1)*thres));
    str=sprintf('%d %d %d',0,length(indu),class); 
    str = [str char(12*ones(1,10-length(str)))];
    fprintf(fid2,'%s\n',str);    
    for u=indu,
        a=num2str(data(u,:),preci_str);
        i=find(~isspace(a));
        a=a(i);
        fprintf(fid2,'%s: %s\n',label{u},a);
    end
    fclose(fid2);
end

fclose(fid);
    
    
    
    