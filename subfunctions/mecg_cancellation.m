function residual = mecg_cancellation(peaks,ecg,method,debug,varargin)
% MECG cancellation algorithms using template subtraction like methods.
% Five template subtraction techniques are implemented going from the least
% adaptive to the more adaptive ones:
% TS,TS-CERUTTI,TS-SUZANNA,TS-LP,TS-PCA. If a more adaptive technique is
% required then an the EKF technique as in Sameni et al. is the best.
%
% inputs
%   peaks:      MQRS markers in ms. Each marker corresponds to the
%               position of a MQRS
%   ecg:        matrix of abdominal ecg channels
%   method:     method to use (TS,TS-CERUTTI,TS-SUZANNA,TS-LP,TS-PCA)
%   varargin:
%       nbCycles:   number of cycles to use in order to build the mean MECG template
%       NbPC:       number of principal components to use for PCA
%       fs:         sampling frequency (NOTE: this code is meant to work at 1kHz)
%
% output
%   residual:   residual containing the FECG
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


% == manage inputs
nbCycles = 20;
NbPC = 2;
fs = 1000;

switch nargin
    case 4
    case 5
        nbCycles = varargin{1};
    case 6
        nbCycles = varargin{1};
        NbPC = varargin{2};        
    case 7
        nbCycles = varargin{1};
        NbPC = varargin{2};        
        fs = varargin{3};
    otherwise
        error('mecg_cancellation: wrong number of input arguments \n');        
end

% check that we have more peaks than nbCycles
if nbCycles>length(peaks)
    error('MECGcancellation Error: more peaks than number of cycles for average ecg');
end

% == constants
NB_CYCLES = nbCycles;
NB_MQRS = length(peaks);
ecg_temp = zeros(1,length(ecg));
ecg_buff = zeros(0.7*fs,NB_CYCLES); % ecg stack buffer
Pstart = 0.25*fs-1;
Tstop = 0.45*fs;

try
    % == template ecg (TECG)
    indMQRSpeaks = find(peaks>Pstart);
    for cc=1:NB_CYCLES
        peak_nb = peaks(indMQRSpeaks(cc+1));   % +1 to unsure full cycles
        ecg_buff(:,cc) = ecg(peak_nb-Pstart:peak_nb+Tstop)';
    end
    TECG = median(ecg_buff,2);

    if strcmp(method,'TS-PCA'); [U,~,~] = svds(ecg_buff,NbPC); end;

    % == MECG cancellation
    for qq=1:NB_MQRS
        if peaks(qq)>Pstart && length(ecg)-peaks(qq)>Tstop

            if strcmp(method,'TS')
                % - simple TS -
                ecg_temp(peaks(qq)-Pstart:peaks(qq)+Tstop) = TECG';
            elseif strcmp(method,'TS-SUZANNA')
                % - three scaling factors - 
                M  = zeros (0.7*fs,3);
                M(1:0.2*fs,1) = TECG(1:Pstart-0.05*fs+1);
                M(0.2*fs+1:0.3*fs,2) = TECG(Pstart-0.05*fs+2:Pstart+0.05*fs+1);
                M(0.3*fs+1:end,3) = TECG(Pstart+2+0.05*fs:Pstart+1+Tstop);
                a = (M'*M)\M'*ecg(peaks(qq)-Pstart:peaks(qq)+Tstop)';
                ecg_temp(peaks(qq)-Pstart:peaks(qq)+Tstop) = a(1)*M(:,1)'+a(2)*M(:,2)'+a(3)*M(:,3)';

            elseif strcmp(method,'TS-CERUTTI')
                % - only one scaling factor -
                M = TECG;
                a = (M'*M)\M'*ecg(peaks(qq)-Pstart:peaks(qq)+Tstop)';
                ecg_temp(peaks(qq)-Pstart:peaks(qq)+Tstop) = a*M';    

            elseif strcmp(method,'TS-LP')
                % - Linear prediction method (Ungureanu et al., 2007) -
                % NOTE: in Ungureanu nbCycles=7
                if qq>NB_CYCLES
                    M = ecg(peaks(qq)-Pstart:peaks(qq)+Tstop)';
                    Lambda = (ecg_buff'*ecg_buff)\ecg_buff'*M;
                    if sum(isnan(Lambda))>0
                        Lambda = ones(length(Lambda),1)/(length(Lambda));
                    end
                    ecg_temp(peaks(qq)-Pstart:peaks(qq)+Tstop) = Lambda'*ecg_buff';
                else
                    M = TECG;
                    ecg_temp(peaks(qq)-Pstart:peaks(qq)+Tstop) = M';
                end            

            elseif strcmp(method,'TS-PCA')
                if mod(qq,10)==0
                    % - to allow some adaptation of the PCA basis -
                    % !!NOTE: this adaption step is slowing down the code!!
                    [U,~,~]   = svds(ecg_buff,NbPC);
                end
                % - PCA method -
                X_out  = ecg(peaks(qq)-Pstart:peaks(qq)+Tstop)*(U*U');
                ecg_temp(peaks(qq)-Pstart:peaks(qq)+Tstop) = X_out;
            end

            if qq>NB_CYCLES
               % adapt template conditional to new cycle being very similar to
               % meanECG to avoid catching artifacts. (not used for PCA method).
               Match = CompareCycles(TECG', ecg(peaks(qq)-Pstart:peaks(qq)+Tstop)',0.8);
               if Match
                   ecg_buff = circshift(ecg_buff,[0 -1]);
                   ecg_buff(:,end) = ecg(peaks(qq)-Pstart:peaks(qq)+Tstop)';
                   TECG = median(ecg_buff,2);
               end
            end       

        % == managing borders
        elseif peaks(qq)<=Pstart
            % - first cycle if not full cycle -
            n = length(ecg_temp(1:peaks(qq)+Tstop)); % length of first pseudo cycle
            m = length(TECG); % length of a pseudo cycle
            ecg_temp(1:peaks(qq)+Tstop+1) = TECG(m-n:end);
        elseif length(ecg)-peaks(qq)<Tstop
            % - last cycle if not full cycle -
            ecg_temp(peaks(qq)-Pstart:end) = TECG(1:length(ecg_temp(peaks(qq)-Pstart:end)));
        end  
    end

    % compute residual
    residual = ecg - ecg_temp;

catch ME
    for enb=1:length(ME.stack); disp(ME.stack(enb)); end;
    residual = ecg;
end

% == debug
if debug
   FONT_SIZE = 15;
   tm = 1/fs:1/fs:length(residual)/fs;
   figure('name','MECG cancellation');
   plot(tm,ecg,'LineWidth',3);
   hold on, plot(tm,ecg-residual,'--k','LineWidth',3);
   hold on, plot(tm,residual-1.5,'--r','LineWidth',3);
   hold on, plot(tm(peaks),ecg(peaks),'+r','LineWidth',2);
   
   legend('mixture','template','residual','MQRS'); 
   title('Template subtraction for extracting the FECG');
   xlabel('Time [sec]'); ylabel('Amplitude [NU]')
   set(gca,'FontSize',FONT_SIZE);
   set(findall(gcf,'type','text'),'fontSize',FONT_SIZE);
end

end

function match = CompareCycles(cycleA,cycleB,thres)
% cross correlation measure to compare if new ecg cycle match with template.
% If not then it is not taken into account for updating the template.
    match = 0;
    bins = size(cycleA,2);
    coeff = sqrt(abs((cycleA-mean(cycleA))*(cycleB-mean(cycleB)))./(bins*std(cycleA)*std(cycleB)));
    if coeff > thres; match = 1; end;
end






















