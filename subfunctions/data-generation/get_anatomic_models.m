function models = get_anatomic_models()
% function models = get_anatomic_models()
% Reads anatomic models currently available in FECGSYN's 
% '/data/anatomic_models/' directory
%
% outputs:
%   models: Struct of available anatomic models
%       models.modelname = '/path/to/model/'
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

conf = FECGSYN_check_platform();

no_models_found_str = ['There are no anatomic models available in your fecgsyn directory.' newline ...
        'Follow the instructions <a href="http://fernandoandreotti.github.io/fecgsyn/pages/install.html#install-models">here</a> to download a set of pre-processed models.'];

% Set anatomic_models directory
anatomic_models_path = [conf.fecgsynpath 'data' conf.slashchar 'anatomic_models' conf.slashchar];

% Print error if no anatomic_models directory exists
if ~exist(anatomic_models_path, 'dir')
    error(no_models_found_str);
end

% Open anatomic models directory to look for subfolders
dirinfo = dir(anatomic_models_path);

% Remove non-directories
dirinfo(~[dirinfo.isdir]) = [];  
dirs = ones(1,length(dirinfo));
for i=1:length(dirinfo)
    if strcmp(dirinfo(i).name,'.') || strcmp(dirinfo(i).name, '..')
    	dirs(i) = 0;
    end
end
dirinfo(~dirs) = [];

% Look for .mat files one sublevel down
mat_files = {};
for K = 1 : length(dirinfo)
  thisdir = dirinfo(K).name;
  found_files = dir([anatomic_models_path conf.slashchar thisdir conf.slashchar '*.mat']);
  if ~isempty(found_files)
      mat_files{1, K} = thisdir;
      mat_files{2, K} = found_files;
  end
end

% Print error if no mat files were found 
if isempty(mat_files)
    error(no_models_found_str);
end

% Check if any mat files contain a fecgsyn anatomic model
% If so, add to the model struct
num_mat_files  = size(mat_files,2);
models = {};
for i=1:num_mat_files
    full_filename = fullfile(mat_files{2,i}.folder, mat_files{2,i}.name);
    name = string(mat_files{2,i}.name);
    name = erase(name,'.mat');
    fieldname = strjoin([ mat_files{1,i} '_' name]);
    fieldname = char(erase(fieldname,' '));
    temp = load(full_filename);
    if (isfield(temp, 'header'))
        header = temp.header;
        if (isfield(header, 'filetype'))
            if strcmp('FECGSYN ANATOMIC MODEL',header.filetype)
                models = setfield(models, fieldname, struct('folder', mat_files{2,i}.folder, 'name', name));
            end
       end
    end
end

% Print error if no files contain a fecgsyn anatomic model
if isempty(models)
    error(no_models_found_str);
end

end