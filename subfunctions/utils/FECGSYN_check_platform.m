function conf = FECGSYN_check_platform()
% function conf = FECGSYN_check_platform()
% Checks the current platform and returns operating configuration
%
% outputs:
%   conf: Struct of operating configuration
%       conf.fecgsynpath = '/path/to/fecgsyn/'
%       conf.slashchar = '/' | '\' 
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

if isunix
    slashchar = '/';
elseif ispc
    slashchar = '\';
else
    error('Platform currently not supported');
end

dbStackOutput = dbstack('-completenames');
pathParts = strsplit(dbStackOutput(1).file, slashchar);

fecgsynpath = '';
for i=1:length(pathParts)-3
    fecgsynpath = strcat(fecgsynpath, pathParts{i}, slashchar);
end

conf = struct('fecgsynpath', fecgsynpath, 'slashchar', slashchar);

end