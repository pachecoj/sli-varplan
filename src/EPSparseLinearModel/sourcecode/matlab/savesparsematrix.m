function dummy=savesparsematrix(spmat,fname)
%SAVESPARSEMATRIX Stores sparse matrix in LHOTSE SimpSparseMatrix format
% SAVESPARSEMATRIX(SPMAT,FNAME) stores sparse matrix into file with
%   name FNAME, using LHOTSE SimpSparseMatrix file format. FNAME
%   can also be a file identifier (as returned by FOPEN) for a file
%   opened in big-endian byte order. In this case, the file is not
%   closed after writing. See 'SimpSparseMatrix::save' in module
%   'matrix' for the file format.

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

if ~issparse(spmat)
  error('SPMAT must be sparse matrix');
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
fwrite(fid,'@SimpSparseMatrix','uchar');
fwrite(fid,0,'int32'); % FF version number 0
fwrite(fid,8,'int32'); % Byte size of double (8)
[m,n]=size(spmat);
[pr,ir,jc]=fst_getsparse(spmat);
fwrite(fid,m,'int32');
fwrite(fid,n,'int32');
fwrite(fid,jc,'int32');
fwrite(fid,ir,'int32');
fwrite(fid,pr,'double');
if ischar(fname)
  fclose(fid);
end
