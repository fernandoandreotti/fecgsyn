% Generate mixture of MECG, FECG and noise
function [mixture,mecg,fecg,noise] = generate_ecg_mixture(debug,SNRfm,SNRmn,varargin)
% generate ecg mixture (mecg, fecg and noise).
%
% inputs
%        debug:     debug [bool]
%        SNRfm:       SNR of fetal signal with respect to maternal signal
%        SNRmn:       SNR of maternal signal compared to background noise
%   <optional>
%   structure as:   <source>.VCG - VCG signal for given source
%                   <source>.H   - Dower-like matrix H (for propagation of dipole) 
%                   <source>.SNR - Gain which is given to source 
%                   Obs: first source is taken as reference for SNR calculus
%                   <source>.type - Maternal (1), Fetal (2) or Noise (3)
% output
%       mixture: mixture of MECG, FECG and noise
%       mecg:    matrix containing projected maternal ECG signal
%       fecg:    cell array containing projected fetal ECG signal(s)
%       noise:   cell array containing projected noise sources
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

% == checking inputs
if nargin<4; error('No source has been given to model generation'); end;

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

% == general
cpt2 = 0; cpt3 = 0; noise = {};
signalf = zeros(NB_FOETUSES*NB_EL,NB_SAMPS);
signaln = zeros(NB_NOISE*NB_EL,NB_SAMPS);

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
            signaln((cpt3-1)*NB_EL+1:cpt3*NB_EL,:) = signal; 
    end
end

% adding maternal signal to mixture
mixture = mecg;

% == SNR calculation for different sources
Pm = sqrt(sum(mixture.^2,2)); % average power of maternal signal throughout channels

% == calibrating FECG (fetal - mother)
% calibration is done using mean maternal and fetal ECG signal powers are reference
if ~isempty(signalf)
    fecg = cell(size(signalf,1)/NB_EL,1);
    powerm = mean(Pm); % mean power of maternal signal accross all channels
    ampf = reshape(sqrt(sum((signalf).^2,2)),NB_EL,[]); % power of each fetus in one column
    powerf = mean(ampf); % mean power of EACH fetal signal (since VCGs are not normalized)
    % run through sources so that every source so that each fetal ECG has SNRfm [dB].
    for i = 1:size(signalf,1)/NB_EL
        % calibrating different hearts
        p = sqrt(exp(-log(10)*-SNRfm/10)*powerm/powerf(i));
        fblock = p*signalf((i-1)*NB_EL+1:i*NB_EL,:);
        mixture = mixture + fblock;
        fecg{i} = fblock;
    end
end

% == calibrating NOISE (noise - fetal)
% re-scale so that together, all noise sources have a 1/SNRmn [dB] level
if ~isempty(signaln)
    noise = cell(size(signaln,1)/NB_EL);
    ampn = reshape(sqrt(sum((signaln).^2,2)),NB_EL,[]); % power of each source in one column (rows are channels)
    ampnorm = diag(1./sum(ampn,2))*ampn; % normalizing in total signal power (%)
    Pnoise = Pm./10^(SNRmn/10); % maximum power allowed per channel
    Peach = diag(Pnoise)*ampnorm; % desired power per channel and per noise source
    % add noise to mixture signals with amplitude modulation
    for i = 1:size(signaln,1)/NB_EL     
        p = 10.^(-SNRmn/20*sqrt(Pm./Peach(:,i)));  % applied gain for each noise signal / channels
        nblock = diag(p)*signaln((i-1)*NB_EL+1:i*NB_EL,:); % re-scaling signals
        mixture = mixture + nblock; % adding noise to mixture
        noise{i} = nblock; % saving signal separetely
    end
end

% == debug
if debug || debug>1 
   FONT_SIZE = 15;
   LINE_WIDTH = 2;
   NB_FET = length(fecg);
   NB_NOISE = length(noise);
   figure
   % plot maternal signal
   subplot(NB_FET+NB_NOISE+2,1,1)
   plot(mecg(1,:),'k','LineWidth',LINE_WIDTH)
   title('Maternal signal ch1')
   set(gca,'FontSize',FONT_SIZE)
   % plot fetal signals
   for i=1:NB_FET
       subplot(NB_FET+NB_NOISE+2,1,1+i)
       plot(fecg{i}(1,:),'b','LineWidth',LINE_WIDTH)
       title(['Fetal signal ' num2str(i) ' ch1'])
       set(gca,'FontSize',FONT_SIZE)
   end
   % plot noises
   for i = 1:NB_NOISE
       subplot(NB_FET+NB_NOISE+2,1,1+NB_FET+i)
       plot(noise{i}(1,:),'r','LineWidth',LINE_WIDTH)
       title(['Noise signal ' num2str(i) ' ch1'])
       set(gca,'FontSize',FONT_SIZE)
   end
   % plot resulting signal
   subplot(NB_FET+NB_NOISE+2,1,NB_FET+NB_NOISE+2)
   plot(mixture(1,:),'m','LineWidth',LINE_WIDTH)
   title('Resulting signal for ch1')
   set(gca,'FontSize',FONT_SIZE);
   xlabel('Sample Number','FontSize',FONT_SIZE)
   set(findall(gcf,'type','text'),'fontSize',FONT_SIZE);   
end

end

