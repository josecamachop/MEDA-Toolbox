function text_tot = cprint(console,text,text_tot,opt,delay,maxl,maxr)

% Print a text on a console.
%
% cprint(console,text) % standard options
% cprint(console,text,text_tot,opt,delay,maxl,maxr) % complete call
%
%
% INPUTS:
%
% console: (1x1) handle of the object, 0 for main console.
%
% text: (text) message to be written.
%
% text_tot: (text) complete text ([] by default).
%
% opt: (1x1) option to refresh the console:
%       - -1: clear console.
%       - 0: clear console and write a new line.
%       - 1: write a new line (by default).
%       - 2: write in the current line.
%
% delay: (1x1) delay in seconds after writting (0.1 by default).
%
% maxl: (1x1) maximum length of a line, in characters (500 by default).
%
% maxr: (1x1) maximum number of lines in console (10 by default).
%
%
% OUTPUTS:
%
% text_tot: (text) output text resulting from the concatenation of the two
%   texts in the input.
%
%
% codified by: José Camacho
% last modification: 08/Jul/15.
%
% Copyright (C) 2016  University of Granada, Granada
% Copyright (C) 2016  José Camacho
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

% Parameters cheking

if nargin < 3, text_tot=' '; end
if nargin < 4, opt=1; end
if nargin < 5, delay=0.1; end
if nargin < 6, maxl=500; end
if nargin < 7, maxr=10; end

% Main code

if console==0,
    disp(text);
    return
end



switch opt,
    case -1,
        set(console,'String',[]);
    case {0,1}
       text = cat(2,sprintf('>>  '),text);
       text_tot = strvcat(text_tot,text);
    case 2,
        text2=cat(2,deblank(text_tot(end,:)),' ',text);
        text_tot = strvcat(text_tot(1:end-1,:),text2);
        textG=get(console,'String');
        text=cat(2,deblank(textG(end,:)),' ',text);
end

rows=floor(length(text)/maxl);
textb=text(1:min(length(text),maxl));
for ind=1:rows,
    textb=strvcat(textb,cat(2,'   ',text(ind*maxl+1:min(length(text),(ind+1)*maxl))));
end
text = textb;

switch opt,
    case 0,
        set(console,'String',text);
    case 1,
        textG=get(console,'String');
        st=size(textG);
        if st(1)>=maxr, textG=textG(2:end,:); end
        set(console,'String',strvcat(textG,text));
    case 2,
        set(console,'String',strvcat(textG(1:end-1,:),text));
end


pause(delay)