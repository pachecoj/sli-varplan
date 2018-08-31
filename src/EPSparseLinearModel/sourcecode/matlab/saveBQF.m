function dummy=saveBQF(fname,umat,a,b,varargin)
%SAVEBQF Stores data matrix in LHOTSE byte-quantized format (BQF)
% SAVEBQF(FNAME,UMAT,A,B,{LABELS},{MAXLVAL}) stores data matrix
%   into BQF file with name FNAME. UMAT is the unsigned byte version
%   of the data matrix, and A,B are the linear conversion factors,
%   i.e. DATA = A*UMAT + B.
%   Optionally, LABELS contains labels which are unsigned byte
%   numbers. In this case, also MAXLVAL must be given and contain
%   the maximum possible label value. All entries of LABELS must
%   lie between 0 and MAXLVAL.
%
%   NOTE: This function cannot produce BQF files with regression
%   targets!

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

fid=fopen(fname,'w');
if fid==-1
  error('Cannot create file!');
end
fprintf(fid,'LHOTSE:BQF\n1\n');
if nargin>4
  if nargin~=6
    error('Require 4 or 6 arguments!');
  end
  writelabels=1;
  labels=varargin{1};
  maxlval=varargin{2};
  if maxlval>255 | maxlval<0
    error('MAXLVAL must be an integer between 0 and 255!');
  end
  if ~isempty(find(labels<0)) | ~isempty(find(labels>maxlval))
    error('LABELS must contain entries between 0 and MAXLVAL!');
  end
else
  writelabels=0;
end
fprintf(fid,'%d\n',writelabels);
if writelabels
  fprintf(fid,'%d\n',maxlval);
end
[numcases numattr]=size(umat);
if writelabels
  [d1 d2]=size(labels);
  if d1~=numcases | d2~=1
    error('LABELS has wrong size!');
  end
end
fprintf(fid,'%d\n',numattr);
fprintf(fid,'%15.14e\n%15.14e\n',a,b);
fprintf(fid,'%d\n',numcases);
fwrite(fid,umat','uchar');
if writelabels
  fwrite(fid,labels,'uchar');
end
fclose(fid);
