function add_data(name,path,data,class,type,thres,preci,debug)

% Add data to a file in the clustering file system. 
%
% add_data(name,path,data,class,type,thres) % minimum call
% add_data(name,path,data,class,type,thres,preci,debug) % complete call
%
%
% INPUTS:
%
% name: (str) name of the file.
%
% path: (str) path to the directory where the clustering data files are
%   located.
%
% data: (LxM) observations to include in the file.
%
% class: (1x1) class associated to the observations.
%
% type: (1xL) type of update:
%       'w'     open file for writing; discard existing contents
%       'a'     open or create file for writing; append data to end of file
%
% thres: (1x1) maximum number of entries in a file.
%
% preci: (1x1) number of decimals (8 by default)
%
% debug: (1x1) disply debug messages
%       0: no messages are displayed.
%       1: display only main messages (default) In the present routine, no 
%           messages are displayed.
%       2: display all messages.
%
%
% OUTPUTS:
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 24/Jan/14.
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
    
%

% Parameters checking 

if nargin < 6, error('Error in the number of arguments.'); end;
if nargin < 7, preci=8; end;
if isempty(preci), preci=8; end;
if nargin < 8, debug = 1; end;

% Computation

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
            data2 = read_data(name,path,s(2),debug);
            data = [data2;data];
            add_data1(name,path,data,class,'w',thres,1,preci);
        else
            fid=fopen(file,'r+');
            str=sprintf('%d %d %d',0,stot,class); 
            str = [str 12*ones(1,10-length(str))];
            fprintf(fid,'%s\n',str);  
            fseek(fid,0,'eof');
            for u=1:s(1),
                a=num2str(data(u,:),preci_str);
                i=find(~isspace(a));
                a=a(i);
                fprintf(fid,'%s\n',a);
            end
            fclose(fid);
        end
    else
        for i=1:s2,
            name2 = fscanf(fid,'%s',1);
        end  
        fclose(fid);
        data2 = read_data(name2,path,s(2),debug);
        data = [data2;data];
        add_data1(name,path,data,class,'a',thres,s2,preci);
    end
else
    if s(1) > thres,
        add_data1(name,path,data,class,'w',thres,1,preci);
    else
        fid=fopen(file,'w');
        str=sprintf('%d %d %d',0,s(1),class); 
        str = [str 12*ones(1,10-length(str))];
        fprintf(fid,'%s\n',str);  
        for u=1:s(1),
            a=num2str(data(u,:),preci_str);
            i=find(~isspace(a));
            a=a(i);
            fprintf(fid,'%s\n',a);
        end
        fclose(fid);
    end
end


function add_data1(name,path,data,class,type,thres,s2,preci)

if nargin < 7, error('Error in the number of arguments.'); end;
if nargin < 8, preci=8; end;

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
str = [str 12*ones(1,10-length(str))];
fprintf(fid,'%s\n',str);
fseek(fid,0,'eof');


for i=s2:nfich+s2-1,
    if ~(i==s2 && isequal('a',type)), fprintf(fid,'%s_%d\n',name,i); end;
    file=[path name '_' num2str(i) '.txt'];
    fid2=fopen(file,'w');
    indu = (((i-s2)*thres+1):min(s,(i-s2+1)*thres));
    str=sprintf('%d %d %d',0,length(indu),class); 
    str = [str 12*ones(1,10-length(str))];
    fprintf(fid2,'%s\n',str);    
    for u=indu,
        a=num2str(data(u,:),preci_str);
        i=find(~isspace(a));
        a=a(i);
        fprintf(fid2,'%s\n',a);
    end
    fclose(fid2);
end

fclose(fid);
    
    
    
    