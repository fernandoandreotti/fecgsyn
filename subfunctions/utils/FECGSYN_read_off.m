function  [vertices, faces] = FECGSYN_read_off(filename)
% [vertices, faces] = FECGSYN_read_off(filename) 
% Reads in the selected .OFF file
%
% inputs:
%   filename: string - filename of input file
%
% outputs:
%   vertices: m x 3 array - vertex coordinates
%   faces: n x 3 array  - vertex connectivity
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
% Last updated : 30-08-2019
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
if ~strcmpi(filetype(1:3), 'off')
    error('This filetype is not supported. Please select an .off file.');    
end

elementInfo = split(fgets(fid));
numVertices = str2double(elementInfo(1));
numFaces = str2double(elementInfo(2));

% Read in vertices
% 3 data points per vertex (x,y,z)
[vertexData,verticeCount] = fscanf(fid,'%f %f %f', 3*numVertices);
if verticeCount~=3*numVertices
     error('File contains incorrect number of vertices.');
end
vertices = reshape(vertexData, 3, verticeCount/3);
vertices = transpose(vertices);

% Read in faces
% 4 data points per face (face number and 3 connected vertices)
[faceData,faceCount] = fscanf(fid,'%d %d %d %d\n', 4*numFaces);
if faceCount~=4*numFaces
     error('File contains incorrect number of faces.');
end
faces = reshape(faceData, 4, faceCount/4);
faces = transpose(faces(2:4,:)+1); % Remove face number and shift

fclose(fid);

end