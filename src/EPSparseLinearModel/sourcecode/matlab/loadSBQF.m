function [data,a,b,varargout]=loadSBQF(fname)
%LOADSBQF Loads data file stored in LHOTSE sparse byte-quantized format (SBQF)
%  [DATA,A,B,{LABELS},{MAXLVAL}] = LOADSBQF(FNAME) loads file with
%  name FNAME. The UCHAR data matrix is returned in DATA, the
%  DOUBLE data matrix can be obtained as A*DOUBLE(DATA)+B. Optionally,
%  the labels vector is returned in LABELS (if the file contains
%  labels). In this case, the max. possible label value can also be
%  returned in MAXLVAL.

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

fid=fopen(fname,'r','ieee-be');
if fid==-1
  error('Cannot open file for reading!');
end
line=fgetl(fid);
if ~strcmp(line,'LHOTSE:SBQF') && ~strcmp(line,'STATSIM:SBQF')
  error('File is not in SBQF format!');
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
  error('Cannot load SBQF files with regression targets!');
end
if nargout<=3
  loadlabels=0;
end
if loadlabels
  maxlval=str2num(fgetl(fid));
end
numattr=str2num(fgetl(fid));
a=str2num(fgetl(fid));
b=str2num(fgetl(fid));
numcases=str2num(fgetl(fid));
%numcases=3; % DEBUG!

data=uint8([]);
for i=1:numcases
  row=zeros(1,numattr);
  num=fread(fid,1,'int');
  if num<0 | num>numattr,
    error('Invalid number of elements for pattern');
  end
  offset=fread(fid,num,'uint8');
  large=find(offset==0); ll=length(large);
  if ll>0,
    lvals=fread(fid,length(large),'int');
    offset(large)=lvals;
  end
  vals=fread(fid,num,'uint8');
  ind=cumsum(offset);
  row(ind)=vals;
  data=[data; uint8(row)];
end

if loadlabels
  varargout(1)={fread(fid,[numcases 1],'uint8')};
  if nargout>4
    varargout(2)={maxlval};
  end
end
