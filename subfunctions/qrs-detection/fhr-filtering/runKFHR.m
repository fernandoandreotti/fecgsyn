function runKFHR(fecg,sqi,qrs_ref,fs)
% run Kalman Filter Heart Rate estimation
% 
% This algorithm is the main routine for estimating HR using the KF
% approach as suggested by Li et al 2008. 
% 
% Input:
%    fecg           Fetal ECG signals (NxM matrix, N = number of channels)
%    sqi            Estimated signal quality for each channel N and segment
%                   S (SxN matrix)
%    qrs_ref        reference QRS detections
%    fs             Sampling rate (in Hz)
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
% (SQI indices)
% Andreotti, F., Gräßer, F., Malberg, H., & Zaunseder, S. (2017). Non-Invasive Fetal ECG Signal Quality Assessment for 
% Multichannel Heart Rate Estimation. IEEE Trans. Biomed. Eng., (in press).
%  
% (FECGSYN Toolbox)
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

% Params
chans = 1:size(sqi,2); % define channels to use
R0 = 1e-3;
Q0 = 1;
p = 1;
TOL = 10; % tolerance for HDR (in bpm)

[fhr,fhr_ref] = prepare4kf(fecg,qrs_ref,fs); % get sqis and references
%% Run Kalman filter
%     fhr_ref - FHR reference
%     fhr     - estimated FHR for each channel
%     sqi     - SQI for each channel
[fhr_final,fhr_estim,innov_estim] = multichan_hr_KF(fhr(:,chans),sqi(:,chans),p,Q0,R0); % multi and single channel results
fhr_weight=sum(fhr(:,chans).*sqi(:,chans)./repmat(sum(sqi(:,chans),2),1,length(chans)),2); % weighted average on SQI values

% removing bad quality segments!!
rem = isnan(fhr_ref)';
if any(rem)
    fhr_ref(rem) = []; fhr_estim(rem,:) = [];  fhr_final(rem,:) = []; fhr(rem,:) = []; fhr_weight(rem) = [];
end

% HDR results
hdr_pre =  1-sum(abs(repmat(fhr_ref,1,length(chans))-fhr(:,chans))>TOL)./length(fhr_ref);
hdr_kfsingle =  1-sum(abs(repmat(fhr_ref,1,length(chans))-fhr_estim)>TOL)./length(fhr_ref);
hdr_weighted = 1-sum(abs(fhr_ref-fhr_weight)>TOL)/length(fhr_ref);
hdr_multikf = 1-sum(abs(fhr_ref-fhr_final)>TOL)/length(fhr_ref);

fprintf('Mean HDR multilead Kalman filter %2.2f\n',hdr_multikf)
fprintf('Mean HDR single lead (%2.2f) and teoretical best result possible (%2.2f)\n',mean(hdr_pre),max(hdr_pre))
fprintf('Mean HDR KF single lead (%2.2f) and teoretical best result possible(%2.2f)\n',mean(hdr_kfsingle),max(hdr_kfsingle))
fprintf('Mean HDR multilead weighted average (%2.2f)\n',hdr_weighted)
fprintf('Mean HDR multilead Kalman filter (%2.2f)\n',hdr_multikf)

% RMSE results

rmse_pre =  sqrt(sum((repmat(fhr_ref,1,length(chans))-fhr(:,chans)).^2)/length(fhr_ref));
rmse_kfsingle =  sqrt(sum((repmat(fhr_ref,1,length(chans))-fhr_estim).^2)/length(fhr_ref));
rmse_weighted = sqrt(sum((fhr_ref-fhr_weight).^2)/length(fhr_ref));
rmse_multikf = sqrt(sum((fhr_ref-fhr_final).^2)/length(fhr_ref));
fprintf('Mean RMSE single lead (%2.2f) and teoretical best result possible (%2.2f)\n',mean(rmse_pre),min(rmse_pre))
fprintf('Mean RMSE KF single lead (%2.2f) and teoretical best result possible (%2.2f)\n',mean(rmse_kfsingle),min(rmse_kfsingle))
fprintf('Mean RMSE multilead weighted average (%2.2f)\n',rmse_weighted)
fprintf('Mean RMSE multilead Kalman filter (%2.2f)\n',rmse_multikf)


% Plots
h = figure;
ax(1)=subplot(3,1,1);
plot(fhr(:,chans),'*')
hold on
set(gca,'ColorOrderIndex',1)
plot(fhr_estim)
plot(fhr_ref,'k','LineWidth',2)
plot(fhr_final,'--g','LineWidth',2)
ylim([100 200])
legend('Rough FHR_1','Rough FHR_2','Rough FHR_3','KF FHR_1','KF FHR_2','KF FHR_3','FHR_{ref}','FHR_{final}')
ylabel('FHR (in bpm)')
ax(2)=subplot(3,1,2);
set(gca,'ColorOrderIndex',1)
plot(100.*sqi(:,chans))
ylim([0 100])
ylabel('SQI (in %)')
ax(3)=subplot(3,1,3);
set(gca,'ColorOrderIndex',1)
plot(innov_estim)
ylabel('Innovation (\nu)')
linkaxes(ax,'x')
%                 xlim([990 1050])
%         saveas(h,sprintfsqinorm('FHR%d_chan%d.png',i,c)); % will create FIG1, FIG2,...

end

function [fhr,fhr_ref] = prepare4kf(fecg,qrs_ref,fs)
% Params
win0 = 5;   % window size for[fhr,fhr_ref] = prepare4kf(fecg,qrs_ref,fs); % get sqis and references
olap0 = 1;  % overalaping [s]
NCHAN = size(fecg,1);
win = win0*fs;
olap = olap0*fs;

%% Obtaining segment-wise FHR, references and
k = 1;
seg = 1;
fhr = zeros(floor(size(fecg,2)/win),NCHAN); % prealloc
fhr_ref = zeros(floor(size(fecg,2)/win),1); % prealloc
fqrs = cell(size(fecg,1),1); % prealloc

for j = 1:size(fecg,1)
    fqrs{j} = OSET_MaxSearch(fecg(j,:),1.3/fs); % 2.3 Hz is the expected FHR parameter
end


while k < size(fecg,2)-win       % Loop through segments
    fref = qrs_ref(qrs_ref>k&qrs_ref<k+win)-k;
    for ch = 1:NCHAN                 % Loop through channels
        fqrsloc = fqrs{ch}(fqrs{ch}>=k&fqrs{ch}<k+win)-k+1;
        fhr(seg,ch) = fs*60/mean(diff(fqrsloc));  % in bpm
    end
    %sqi(seg,:) = table2array(Data.sqi(idx,11)); % bsqi13
    fhr_ref(seg,1) =  fs*60/mean(diff(fref));
    k = k+olap;  % move window
    seg = seg+1; % loop table indexes
    %disp(['Sample: ' num2str(k)])
end
% using sigmoid function to fit SQI "classes" into regression curve

end
