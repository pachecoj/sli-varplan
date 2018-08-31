function dummy=savebasematrix(mat,fname,varargin)
%SAVEBASEMATRIX Stores matrix in LHOTSE BaseMatrix format
% SAVEBASEMATRIX(MAT,FNAME,{TRUETYPE=6},{NOTAG=0}) stores matrix
%   MAT into file with name FNAME, using LHOTSE BaseMatrix file
%   format. FNAME can also be a file identifier (as returned by
%   FOPEN) for a file opened in big-endian byte order. In this
%   case, the file is not closed after writing. If TRUETYPE is
%   given, MAT is converted into the following type before writing:
%   - 1: schar (signed char)
%   - 2: uchar (unsigned char)
%   - 3: int32
%   - 4: uint32
%   - 5: float32
%   - 6: double [default]
%   - 7: uint8 (this is bool in LHOTSE)
%   If NOTAG~=0, we do not write a file tag.
%
%   NOTE: The format of LHOTSE StMatrix is the same as
%   BaseMatrix<double>.

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

truetype=6; notag=0;
if nargin>2,
  truetype=varargin{1};
  if nargin>3,
    notag=(varargin{2}~=0);
  end
end
if truetype<1 | truetype>7
  error('Unsupported TRUETYPE type code');
end
[m,n]=size(mat);
% Create file in big-endian byte order (LHOTSE default)
if ischar(fname)
  fid=fopen(fname,'w','ieee-be');
  if fid==-1
    error(['Cannot create ''' fname '''!']);
  end
else
  fid=fname;
end
if ~notag,
  fwrite(fid,'@BaseMatrix','uchar');
end
fwrite(fid,1,'int32'); % FF version number 1
fwrite(fid,truetype,'int32');
fwrite(fid,m,'int32');
fwrite(fid,n,'int32');
if truetype==6,
  fwrite(fid,8,'int32');
  fwrite(fid,mat,'double');
elseif truetype==1,
  fwrite(fid,1,'int32');
  fwrite(fid,mat,'schar');
elseif truetype==2,
  fwrite(fid,1,'int32');
  fwrite(fid,mat,'uchar');
elseif truetype==3,
  fwrite(fid,4,'int32');
  fwrite(fid,mat,'int32');
elseif truetype==4,
  fwrite(fid,4,'int32');
  fwrite(fid,mata,'uint32');
elseif truetype==5,
  fwrite(fid,4,'int32');
  fwrite(fid,mat,'float32');
else
  fwrite(fid,1,'int32');
  fwrite(fid,mat,'uint8');
end
if ischar(fname)
  fclose(fid);
end
