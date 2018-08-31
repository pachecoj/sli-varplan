function mat=loadbasematrix(fname,varargin)
%LOADBASEMATRIX Loads matrix stored in LHOTSE BaseMatrix format
% MAT=LOADBASEMATRIX(FNAME,{TRUETYPE=-1},{RET=0},{NOTAG=0}) loads matrix
%   from file with name FNAME, written in LHOTSE BaseMatrix file format.
%   The matrix is returned in MAT. The file format contains type
%   information. FNAME can also be a file identifier (as returned by
%   FOPEN) for a file opened in big-endian byte order. In this case, the
%   file is not closed after reading.
%   If TRUETYPE is given, the file type code is compared against TRUETYPE,
%   and an exception is thrown if they are different. If TRUETYPE==-1, it
%   is ignored.
%   By default, MAT is a DOUBLE matrix, but if RET is given and ~=0, MAT
%   has the type of the file.
%   If NOTAG~=0, we do not read the file tag.
%
%   NOTE: LHOTSE StMatrix uses BaseMatrix<double> format.
%   DEBUG: Type codes of matrices written by LHOTSE are invalid. If
%   TRUETYPE is given, we use this type to load the matrix, ignoring the
%   typecode field in the file.

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

truetype=-1; ret=0; notag=0;
if nargin>1,
  truetype=varargin{1};
  if nargin>2,
    ret=(varargin{2}~=0);
    if nargin>3,
      notag=(varargin{3}~=0);
    end
  end
end
% Open file in big-endian byte order (LHOTSE default)
if ischar(fname)
  fid=fopen(fname,'r','ieee-be');
  if fid==-1
    error(['Cannot open ''' fname ''' for reading!']);
  end
else
  fid=fname;
end
if ~notag,
  tag=char(fread(fid,[1 11],'uchar'));
  if ~strcmp(tag,'@BaseMatrix'),
    error('Invalid tag!');
  end
end
ffver=fread(fid,1,'int32');
if ffver~=0 && ffver~=1
  error('Unsupported file format number!');
end
tc=fread(fid,1,'int32');
if truetype==-1
  if tc<=0 | tc>7,
    error('Unsupported type code!');
  end
else
  tc=truetype; % DEBUG: Do not use typecode field
end
if tc==1,
  tcs='schar';
elseif tc==2,
  tcs='uchar';
elseif tc==3,
  tcs='int32';
elseif tc==4,
  tcs='uint32';
elseif tc==5,
  tcs='float32';
elseif tc==6,
  tcs='double';
elseif tc==7,
  error('Type bool not supported in the moment!');
end
m=fread(fid,1,'int32');
n=fread(fid,1,'int32');
if n==0 | m==0
  mat=[];
else
  if ffver==1
    % File byte size. Can only deal with fixed sizes
    bsize=fread(fid,1,'int32');
    if (tc==1 && bsize~=1) || (tc==2 && bsize~=1) || ...
	   (tc==3 && bsize~=4) || (tc==4 && bsize~=4) || ...
	   (tc==5 && bsize~=4) || (tc==6 && bsize~=8)
      error('Cannot read file with non-standard byte size');
    end
  end
  if ret,
    tcs=strcat(tcs,'=>',tcs);
  end
  mat=reshape(fread(fid,m*n,tcs),m,n);
end
if ischar(fname)
  fclose(fid);
end
