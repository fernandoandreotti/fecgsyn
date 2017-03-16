function [mixture,mecg,fecg,noise,f_handle] = ...
    generate_ecg_mixture(debug,SNRfm,SNRmn,mqrs,fqrs,fs,varargin)
% function [mixture,mecg,fecg,noise,f_handle] = generate_ecg_mixture(debug,SNRfm,SNRmn,mqrs,fqrs,fs,varargin)
% generate ecg mixture (mecg, fecg and noise).
%
% Input:
%        debug:      debug [bool]
%        SNRfm:      SNR of fetal signal with respect to maternal signal
%        SNRmn:      SNR of maternal signal compared to background noise
%        mqrs:       maternal qrs locations
%        fqrs:       foetal qrs locations
%   <optional>
%   structure as:   <source>.VCG - VCG signal for given source
%                   <source>.H   - Dower-like matrix H (for propagation of dipole) 
%                   <source>.SNR - Gain which is given to source 
%                   Obs: first source is taken as reference for SNR calculus
%                   <source>.type - Maternal (1), Fetal (2) or Noise (3)
% Output:
%        mixture: mixture of MECG, FECG and noise
%        mecg:    matrix containing projected maternal ECG signal
%        fecg:    cell array containing projected fetal ECG signal(s)
%        noise:   cell array containing projected noise sources
%
%
% Examples:
% TODO
%
% See also:
% run_ecg_generator
% 
% 
% --
% fecgsyn toolbox, version 1.2, March 2017
% Released under the GNU General Public License
%
% Copyright (C) 2017  Joachim Behar & Fernando Andreotti
% Department of Engineering Science, University of Oxford
% joachim.behar@oxfordalumni.org, fernando.andreotti@eng.ox.ac.uk
%
% 
% For more information visit: http://www.fecgsyn.com
% 
% Referencing this work
%
% Behar, J., Andreotti, F., Zaunseder, S., Li, Q., Oster, J., & Clifford, G. D. (2014). An ECG Model for Simulating 
% Maternal-Foetal Activity Mixtures on Abdominal ECG Recordings. Physiol. Meas., 35(8), 1537–1550.
% 
%
% Last updated : 15-03-2017
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

% == checking inputs
if nargin<5; error('No source has been given to model generation'); end;
f_handle = [];
mecg = [];
fecg = [];
noise = {};

% == constants
NB_EL = size(varargin{1}.H,1); % number of electrodes
NB_SIG2MIX = length(varargin); % number of signals to mix
NB_SAMPS = size(varargin{1}.VCG,2); % number of samples = signal length
NB_FOETUSES = 0;
NB_NOISE = 0;
for vv=1:NB_SIG2MIX;  
    if (varargin{vv}.type==2); NB_FOETUSES = NB_FOETUSES+1; end; 
    if (varargin{vv}.type==3); NB_NOISE = NB_NOISE+1; end;
end;
% constants to  help normalize the calibration procedure regarding the
% number of beats present on the signal.
MHR = 60; %     [in bpm]
FHR = 120; %    [in bpm]

% == general
cpt2 = 0; cpt3 = 0;
signalf = zeros(NB_FOETUSES*NB_EL,NB_SAMPS);
signaln = {};

% == loop through input dipole structures
for i=1:NB_SIG2MIX
    src = varargin{i}; % src = source structure for input source
    if ndims(src.H)==3 % for time varying H
        VCG = num2cell(src.VCG,1); % transforming VCG matrix to use cellfun
        H = num2cell(src.H,[1 2]); % transforming H matrix to use cellfun
        H = reshape(H,size(VCG)); % resizing matrix (loosing one dimension to match VCG)
        signal = cellfun(@(A,B) A*B,H,VCG,'UniformOutput',false); % projecting sources H(t)*VCG
        signal = cell2mat(signal); % converting back to matrix array
    else % for time invariant H
        signal = src.H*src.VCG;
    end
    % resulting projected signals for ith source
    switch src.type
        case 1 % maternal source
            mecg = signal;          
        case 2 % fetal source
            cpt2 = cpt2+1;
            signalf((cpt2-1)*NB_EL+1:cpt2*NB_EL,:) = signal; 
        case 3 % noise source
            cpt3 = cpt3+1; 
            signal = bsxfun(@times,signal,src.SNRfct); % multiplying by modulating function
            signaln{end+1} = signal; 
    end
end

% adding maternal signal to mixture
mixture = mecg;

% == SNR calculation for different sources
mbeats = 60*fs*length(mqrs)/length(mixture); % now im bpm
Pm = sum(mixture.^2,2)*(MHR/mbeats); % average power of maternal 
                                     % across channels with heart rate correction
powerm = mean(Pm);

% == calibrating FECG (fetal - mother)
% calibration is done using mean maternal and fetal ECG signal powers are reference
if ~isempty(signalf)
    fecg = cell(size(signalf,1)/NB_EL,1);
    fbeats = cellfun(@(x) length(x),fqrs); % multiple foetuses support
    fbeats = 60*fs*fbeats/length(signalf);
    ampf = reshape(sum((signalf).^2,2),NB_EL,[])*diag(FHR./fbeats); % power of each fetus in one column
    powerf = mean(ampf); % mean power of EACH fetal signal 
                         % (since VCGs are not normalized)
     % run through sources so that every source so that each fetal ECG has SNRfm [dB].
    for i = 1:size(signalf,1)/NB_EL
        % calibrating different hearts
        p = sqrt(powerm./powerf(i))*10.^(SNRfm/20);
        fblock = p*signalf((i-1)*NB_EL+1:i*NB_EL,:);
        mixture = mixture + fblock;
        fecg{i} = fblock;
    end
end

% == calibrating NOISE (maternal - noise)
% re-scale so that together, all noise sources have a 1/SNRmn [dB] level
if ~isempty(signaln)
    noise = cell(size(signaln)); % preallocating
    noisegain = cellfun(@(x) mean(sum(x.^2,2)),signaln); % for noises with different power
    noisegain = noisegain./sum(noisegain); % percentual power of each noise source
    
    sig = cat(3,signaln{:});            % transforms cell in 3D matrix
    sigpow = sum(sum(sig,3).^2,2);      % total noise power for each channel
    meannoisepow = mean(sigpow);        % average noise power
    
    p = sqrt(powerm./meannoisepow)*10.^(-SNRmn/20); % applied gain for each noise signal / channels

    % add noise to mixture signals with amplitude modulation
    for i = 1:length(signaln)
        nblock = noisegain(i)*diag(p)*signaln{i}; % re-scaling signals
        mixture = mixture + nblock; % adding noise to mixture
        noise{i} = nblock; % saving signal separetely
    end
end

% == debug
if debug || debug>1 
   FONT_SIZE = 15;
   FONT_SIZE_SMALL = 9;
   LINE_WIDTH = 2;
   NB_FET = length(fecg);
   NB_NOISE = length(noise);
   
   if debug ~= 11
       f_handle = figure('name', 'Generated ECG mixture');
       % plot maternal signal
       subplot(NB_FET+NB_NOISE+2,1,1)
       plot(mecg(1,:),'k','LineWidth',LINE_WIDTH)
       %%%title('Maternal signal ch1')
       ylabel('Maternal signal ch1')
   else
       f_handle(end+1) = figure('name', 'Maternal signal ch1');
       set(f_handle(end), 'Visible', 'off')
       plot(mecg(1,:),'k','LineWidth',LINE_WIDTH)
       xlabel('Sample Number','FontSize',FONT_SIZE_SMALL)
   end
   
   
   %%%set(gca,'FontSize',FONT_SIZE)
   % plot fetal signals
   for i=1:NB_FET
       if debug ~= 11
           subplot(NB_FET+NB_NOISE+2,1,1+i)
           plot(fecg{i}(1,:),'b','LineWidth',LINE_WIDTH)
           %%%title(['Fetal signal ' num2str(i) ' ch1'])
           ylabel(['Fetal signal ' num2str(i) ' ch1'])
           %%%set(gca,'FontSize',FONT_SIZE)
       else
           f_handle(end+1) = figure('name', ['Fetal signal ' num2str(i) ' ch1']);
           set(f_handle(end), 'Visible', 'off')
           plot(fecg{i}(1,:),'b','LineWidth',LINE_WIDTH)
           xlabel('Sample Number','FontSize',FONT_SIZE_SMALL)
       end
   end
   % plot noises
   for i = 1:NB_NOISE
       if debug ~= 11
           subplot(NB_FET+NB_NOISE+2,1,1+NB_FET+i)
           plot(noise{i}(1,:),'r','LineWidth',LINE_WIDTH)
           %%%title(['Noise signal ' num2str(i) ' ch1'])
           ylabel(['Noise signal ' num2str(i) ' ch1'])
           %%%set(gca,'FontSize',FONT_SIZE)
       else
           f_handle(end+1) = figure('name', ['Noise signal ' num2str(i) ' ch1']);
           set(f_handle(end), 'Visible', 'off')
           plot(noise{i}(1,:),'r','LineWidth',LINE_WIDTH)
           xlabel('Sample Number','FontSize',FONT_SIZE_SMALL)
       end
   end
   
   % plot resulting signal
   if debug ~= 11
       subplot(NB_FET+NB_NOISE+2,1,NB_FET+NB_NOISE+2)
       plot(mixture(1,:),'m','LineWidth',LINE_WIDTH)
       %%%title('Resulting signal for ch1')
       ylabel('Resulting signal for ch1')
       %%%set(gca,'FontSize',FONT_SIZE);
       xlabel('Sample Number','FontSize',FONT_SIZE)
       %%%set(findall(gcf,'type','text'),'fontSize',FONT_SIZE);  
   else
       f_handle(end+1) = figure('name', 'Resulting signal for ch1');
       set(f_handle(end), 'Visible', 'off')
       plot(mixture(1,:),'m','LineWidth',LINE_WIDTH)
       xlabel('Sample Number','FontSize',FONT_SIZE_SMALL)
       
   end
end

end

