function [maternalSource, fetalSource] = FECGSYN_read_sources(filename) 
% [maternalSource, fetalSource] = FECGSYN_read_sources(filename) 
% Reads in the selected .sources file
%
% inputs:
%   filename: string - filename of input file
%
% outputs:
%   maternalSource: 1 x 3 array - maternal source coordinates
%   fetalSource: 1 x 3 array - fetal source coordinates
%
% --
% fecgsyn toolbox, version 1.3-alpha, August 2019
% Released under the GNU General Public License
%
% Copyright (C) 2019  Emerson Keenan
% Department of Electrical and Electronic Engineering, University of
% Melbourne
% e.keenan@ieee.org
% 
% For more information visit: http://www.fecgsyn.com
% 
% Referencing this work
%
% Keenan E., Karmakar C K. and Palaniswami M., The effects of asymmetric volume conductor modeling on non-invasive fetal ECG extraction. 
% Physiol Meas 39(10), pp. 105013, 2018.
% 
%
% Last updated : 05-09-2019
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

fid = fopen(filename,'r');

if (fid == -1)
    error('Unable to open the selected file.');
end

filetype = fgets(fid); 
if ~strcmpi(filetype(1:7), 'sources')
    error('This filetype is not supported. Please select a .sources file.');    
end

% Read in maternal source data
% 1 char identifier + x,y,z coordinates
[maternalSourceData,readCount] = fscanf(fid,'%c %f %f %f\n', 4);
if readCount~=4 || maternalSourceData(1) ~= 'M'
     error('Maternal source contains invalid data.');
end
maternalSource = reshape(maternalSourceData(2:4),1,3);

% Read in fetal source data
% 1 char identifier + x,y,z coordinates
[fetalSourceData,readCount] = fscanf(fid,'%c %f %f %f', 4);
if readCount~=4 || fetalSourceData(1) ~= 'F'
     error('Fetal source contains invalid data.');
end
fetalSource = reshape(fetalSourceData(2:4),1,3);

fclose(fid);

end