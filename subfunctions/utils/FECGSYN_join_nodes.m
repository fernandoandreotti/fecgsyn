function [vertices,faces] = FECGSYN_join_nodes(vert1,face1,vert2,face2)
% [vertices,faces] = FECGSYN_join_nodes(vert1,face1,vert2,face2)
% This function takes two sets of input triangular vertices and
% faces, merges them, reindexes duplicate nodes and outputs a set 
% of combined vertices/faces. It assumed that the inputs have 
% no duplicate nodes or intersecting faces within themselves 
% and no intersecting faces with each other
%
% inputs:
%   vert1: m x 3 array - vertices of first mesh
%   face1: n x 3 array - faces of first mesh
%   vert2: m x 3 array - vertices of second mesh
%   face2: n x 3 array - faces of second mesh
%
% outputs:
%   vertices: o x 3 array - vertices of joined mesh
%   faces: o x 3 array - faces of joined mesh
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
% Last updated 05-09-2019
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

if (size(vert1,2) ~= 3 || size(vert2,2) ~= 3)
    error('Input vertices must have 3 columns each');
end
if (size(face1,2) ~= 3 || size(face2,2) ~= 3)
    error('Input faces must have 3 columns each');
end

nVert1 = size(vert1,1);
nVert2 = size(vert2,1);
nFace2 = size(face2,1);
vertices = vert1;
faces = face1;
indexMap = zeros(nVert2,1);

for i=1:nVert2
    vertIndex = -1;
    for j=1:nVert1
       if vert2(i,:) == vert1(j,:)
           vertIndex = j;
       end
    end
    if (vertIndex > 0)
        indexMap(i) = vertIndex;
    else
        vertices = [vertices; vert2(i,:)];
        vertIndex = size(vertices,1);
        indexMap(i) = vertIndex;
    end
end

for i=1:nFace2
    newFaces = zeros(1,3);
    for j=1:3
        newFaces(1,j) = indexMap(face2(i,j));
    end
    faces = [faces; newFaces];
end

end