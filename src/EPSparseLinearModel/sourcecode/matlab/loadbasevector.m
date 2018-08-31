function vec=loadbasevector(fname,varargin)
%LOADBASEVECTOR Loads vector stored in LHOTSE BaseVector format
% VEC=LOADBASEVECTOR(FNAME,{TRUETYPE=-1},{RET=0},{NOTAG=0},{DCTT=-1})
%   loads vector from file with name FNAME, written in LHOTSE
%   BaseVector file format.
%   The (column) vector is returned in VEC. The file format contains type
%   information. FNAME can also be a file identifier (as returned by
%   FOPEN) for a file opened in big-endian byte order. In this case, the
%   file is not closed after reading.
%   If TRUETYPE is given, the file type code is compared against TRUETYPE,
%   and an exception is thrown if they are different. If TRUETYPE==-1, it
%   is ignored. See SAVEBASEVECTOR for type code values.
%   By default, VEC is a DOUBLE vector, but if RET is given and ~=0, VEC
%   has the type of the file.
%   If NOTAG~=0, we do not read the file tag.
%
%   NOTE: LHOTSE StVector uses BaseVector<double> format,
%   IndexVector uses BaseVector<int> format.
%
%   Downward compatibility: Files written by earlier LHOTSE
%   versions come with wrong type code field. If DCTT is given, the
%   file is read using this type code, ignoring the type code field
%   and TRUETYPE.

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

truetype=-1; ret=0; notag=0; dctt=-1;
if nargin>1,
  truetype=varargin{1};
  if nargin>2,
    ret=(varargin{2}~=0);
    if nargin>3,
      notag=(varargin{3}~=0);
      if nargin>4,
	dctt=varargin{4};
      end
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
  %fprintf(1,'%s\n',tag);
  if ~strcmp(tag,'@BaseVector'),
    error('Invalid tag!');
  end
end
ffver=fread(fid,1,'int32');
if ffver~=0 && ffver~=1
  error('Unsupported file format number!');
end
tc=fread(fid,1,'int32');
if dctt~=-1
  tc=dctt;
elseif truetype~=-1 && tc~=truetype
  error('Type code in file different from TRUETYPE');
end
if tc==1,
  tcs='schar'; bs=1;
elseif tc==2,
  tcs='uchar'; bs=1;
elseif tc==3,
  tcs='int32'; bs=4;
elseif tc==4,
  tcs='uint32'; bs=4;
elseif tc==5,
  tcs='float32'; bs=4;
elseif tc==6,
  tcs='double'; bs=8;
elseif tc==7,
  tcs='uint8'; bs=1;
else
  error(sprintf('Unsupported type code %d!',tc));
end
n=fread(fid,1,'int32');
if n==0
  vec=[];
else
  if ffver==1
    % File byte size. Can only deal with fixed sizes
    bsize=fread(fid,1,'int32');
    if bsize~=bs
      error('Cannot read file with non-standard byte size');
    end
  end
  if ret,
    tcs=strcat(tcs,'=>',tcs);
  end
  %tcs
  vec=fread(fid,n,tcs);
end
if ischar(fname)
  fclose(fid);
end
