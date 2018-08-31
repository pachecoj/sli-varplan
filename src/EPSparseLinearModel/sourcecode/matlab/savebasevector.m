function dummy=savebasevector(vec,fname,varargin)
%SAVEBASEVECTOR Stores vector in LHOTSE BaseVector format
% SAVEBASEVECTOR(VEC,FNAME,{TRUETYPE=6},{NOTAG=0}) stores vector
%   VEC into file with name FNAME, using LHOTSE BaseVector file format.
%   FNAME can also be a file identifier (as returned by FOPEN) for a
%   file opened in big-endian byte order. In this case, the file is
%   not closed after writing. If TRUETYPE is given, VEC is converted
%   into the following type before writing:
%   - 1: schar (signed char)
%   - 2: uchar (unsigned char)
%   - 3: int32
%   - 4: uint32
%   - 5: float32
%   - 6: double [default]
%   - 7: uint8 (this is bool in LHOTSE)
%   If NOTAG~=0, we do not write a file tag.
%
%   NOTE: The format of LHOTSE StVector is the same as
%   BaseVector<double>.

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
[n,m]=size(vec);
if n==1,
  n=m;
elseif m~=1,
  error('VEC is not a vector');
end
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
  fwrite(fid,'@BaseVector','uchar');
end
fwrite(fid,1,'int32'); % FF version number 1
fwrite(fid,truetype,'int32');
fwrite(fid,n,'int32');
if truetype==6,
  fwrite(fid,8,'int32');
  fwrite(fid,vec,'double');
elseif truetype==1,
  fwrite(fid,1,'int32');
  fwrite(fid,vec,'schar');
elseif truetype==2,
  fwrite(fid,1,'int32');
  fwrite(fid,vec,'uchar');
elseif truetype==3,
  fwrite(fid,4,'int32');
  fwrite(fid,vec,'int32');
elseif truetype==4,
  fwrite(fid,4,'int32');
  fwrite(fid,vec,'uint32');
elseif truetype==5
  fwrite(fid,4,'int32');
  fwrite(fid,vec,'float32');
else
  fwrite(fid,1,'int32');
  fwrite(fid,vec,'uint8');
end
if ischar(fname)
  fclose(fid);
end
