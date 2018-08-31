function [data,varargout]=loadBQF(fname)
%LOADBQF Loads data file stored in LHOTSE byte-quantized format (BQF)
% [DATA,{LABELS},{MAXLVAL}] = LOADBQF(FNAME) loads file with name FNAME. The
%   data matrix is returned in DATA. Optionally, the labels vector
%   is returned LABELS (if the file contains labels). In this case,
%   the max. possible label value can also be returned in MAXLVAL.

% Copyright (C) 1999-2006 Matthias Seeger
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
% USA.

fid=fopen(fname,'r');
if fid==-1
  error('Cannot open file for reading!');
end
line=fgetl(fid);
if ~strcmp(line,'LHOTSE:BQF') && ~strcmp(line,'STATSIM:BQF')
  error('File is not in BQF format!');
end
ffver=str2num(fgetl(fid));
if ffver~=1
  error('Unknown file format version number!');
end
loadlabels=str2num(fgetl(fid));
if loadlabels<0 | loadlabels>2
  error('Unknown target status!');
end
if loadlabels==2
  error('Cannot load BQF files with regression targets!');
end
if nargout<=1
  loadlabels=0;
end
if loadlabels
  maxlval=str2num(fgetl(fid));
end
numattr=str2num(fgetl(fid));
a=str2num(fgetl(fid));
b=str2num(fgetl(fid));
numcases=str2num(fgetl(fid));
data=a*fread(fid,[numattr numcases],'uchar')'+b;
if loadlabels
  varargout(1)={fread(fid,[numcases 1],'uchar')};
  if nargout>2
    varargout(2)={maxlval};
  end
end
