%% Main script for FECG morphological analysis
% 
% This is the mains script for testing the morphological consistency of
% extracting the foetal signal using various methods. 
% Used extraction methods:
%  - ICA
%  - PCA
%  -piCA
% 
% Used morphological measures:
% - T/QRS ratio
% - ST segment
% - QT interval
% 
% 
% 
% NI-FECG simulator toolbox, version 1.0, February 2014
% Released under the GNU General Public License
%
% Copyright (C) 2014  Joachim Behar & Fernando Andreotti
% Oxford university, Intelligent Patient Monitoring Group - Oxford 2014
% joachim.behar@eng.ox.ac.uk, fernando.andreotti@mailbox.tu-dresden.de
%
% Last updated : 03-06-2014
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

%% Input parameters
% saving path
if isunix
%     path = ['/media/fernando/Data/foobar/fecgdata_test/' datestr(date,'yyyy.mm.dd') '/'];
 path = '/media/fernando/FetalEKG/2014.06_fecgsyn_simulations';
else
    path = ['C:\foobar\fecgdata\' datestr(date,'yyyy.mm.dd') '\'];
end

%% Set-up parameters
generate = 0;   % boolean, data should be generated? 
                % If not, path should direct to data location
                
ch = 1:32;      % channels to be used in ICA

%% Data Generation
if generate
    mkdir(path)
    generate_data(path)  % generates set of unique simulated data for testing
else
    cd(path)
end

%% Extraction Methods
fls = dir('*.mat');     % looking for .mat (creating index)
fls =  arrayfun(@(x)x.name,fls,'UniformOutput',false);
for i = 1:length(fls)       
    disp(['Extracting file ' fls{i} '..'])
    % = loading data
    load(fls{i})
    noise = sum(cat(3,out.noise{:}),3);
    if isempty(noise)
        noise = zeros(size(out.mecg));
    end
    mixture = double(out.mecg) + sum(cat(3,out.fecg{:}),3) ...
        + noise;     % re-creating abdominal mixture
    
   
    % = preprocessing channels
    fs = out.param.fs;
    HF_CUT = 100; % high cut frequency
    LF_CUT = 0.7; % low cut frequency
    wo = 60/(fs/2); bw = wo/35;
    [b_lp,a_lp] = butter(5,HF_CUT/(fs/2),'low');
    [b_bas,a_bas] = butter(3,LF_CUT/(fs/2),'high');
    for j=ch
        lpmix = filtfilt(b_lp,a_lp,mixture(j,:));
        mixture(j,:) = filtfilt(b_bas,a_bas,lpmix);
    end
    
    % = using ICA
    
    loopsec = 60;   % in seconds
    [F1,RMS] = ica_extraction(mixture,fs,ch,out.fqrs{1},loopsec);     % extract using ICA
    
    % = using TSc
    
end
%% Morphological Analysis


