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

% channels to be used in ICA                
% ch = 1:32;      
ch = [1:2:8 10:2:16 17:2:24 26:2:32];
debug = 1;

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
stats_ica = zeros(length(fls),4);
stats_tsc = zeros(length(fls),4);
for i = 1:length(fls)   
    disp(['Extracting file ' fls{i} '..'])
    % = loading data
    load(fls{i})
    noise = sum(cat(3,out.noise{:}),3);
    if isempty(noise)
        noise = zeros(size(out.mecg));
    end
    fs = out.param.fs;
    INTERV = round(0.05*fs);    % BxB acceptance interval
    TH = 0.3;                   % detector threshold
    REFRAC = round(.15*fs)/1000; % detector refractory period
    mixture = double(out.mecg) + sum(cat(3,out.fecg{:}),3) ...
        + noise;     % re-creating abdominal mixture
       
    
    %% Experiment 1
    % = preprocessing channels
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
    disp('ICA extraction ..')
    loopsec = 60;   % in seconds
    icasig = ica_extraction(mixture,fs,ch,out.fqrs{1},loopsec);     % extract using IC
    
    % Calculate quality measures
    qrsica = qrs_detect(icasig,TH,REFRAC,fs);
    [F1,RMS,PPV,SE] = Bxb_compare(out.fqrs{1},qrsica,INTERV);
    stats_ica(end+1,:) = [F1,RMS,PPV,SE];
    
    % = using TSc
    disp('TS extraction ..')
    % look for channel with largest SNRfm
    amps = sum(double(out.fecg{1}).^2,2)./sum(double(out.mecg).^2,2);
    [~,chts]=max(amps);       % chosing the channel with highest fetal signal ratio
    residual = mecg_cancellation(out.mqrs,mixture(chts,:),'TS-CERUTTI',debug);
    qrsts = qrs_detect(residual,TH,REFRAC,fs);
    [F1,RMS,PPV,SE] = Bxb_compare(out.fqrs{1},qrsts,INTERV);
    stats_tsc(end+1,:) = [F1,RMS,PPV,SE];

    % Debug plots
    if debug
        hold on
        plot(out.fqrs{1}/fs,2000,'og','MarkerSize',7)
    end
    clearvars -except stats_ica stats_tsc fls ch debug

end
%% Morphological Analysis


