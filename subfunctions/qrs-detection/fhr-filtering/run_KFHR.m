function run_KFHR(qrs,sqi,qrs_ref)

% Params
chans = 1:size(sqi,2); % define channels to use
R0 = 1e-3;
Q0 = 1;
p = 1;


    [fhr,fhr_ref] = prepare4kf(qrs,qrs_ref); % get sqis and references
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
    hdr_pre =  1-sum(abs(repmat(fhr_ref,1,length(chans))-fhr(:,chans))>tol)./length(fhr_ref);
    hdr_kfsingle =  1-sum(abs(repmat(fhr_ref,1,length(chans))-fhr_estim)>tol)./length(fhr_ref);
    hdr_weighted = 1-sum(abs(fhr_ref-fhr_weight)>tol)/length(fhr_ref);
    hdr_multikf = 1-sum(abs(fhr_ref-fhr_final)>tol)/length(fhr_ref);
    
        sprintf('Mean HDR multilead Kalman filter %2.2f',hdr_multikf)

    sprintf('Mean HDR single lead (%2.2f) and teoretical best result possible (%2.2f)\n',mean(hdr_pre),max(hdr_pre))
    sprintf('Mean HDR KF single lead (%2.2f) and teoretical best result possible(%2.2f)\n',mean(hdr_kfsingle),max(hdr_kfsingle))
    sprintf('Mean HDR multilead weighted average (%2.2f)\n',hdr_weighted)
    sprintf('Mean HDR multilead Kalman filter (%2.2f)\n',hdr_multikf)

    % RMSE results
    
    rmse_pre =  sqrt(sum((repmat(fhr_ref,1,length(chans))-fhr(:,chans)).^2)/length(fhr_ref));
    rmse_kfsingle =  sqrt(sum((repmat(fhr_ref,1,length(chans))-fhr_estim).^2)/length(fhr_ref));
    rmse_weighted = sqrt(sum((fhr_ref-fhr_weight).^2)/length(fhr_ref));
    rmse_multikf = sqrt(sum((fhr_ref-fhr_final).^2)/length(fhr_ref));
    sprintf('Mean RMSE single lead (%2.2f) and teoretical best result possible (%2.2f)\n',mean(rmse_pre),min(rmse_pre))
    sprintf('Mean RMSE KF single lead (%2.2f) and teoretical best result possible (%2.2f)\n',mean(rmse_kfsingle),min(rmse_kfsingle))
    sprintf('Mean RMSE multilead weighted average (%2.2f)\n',rmse_weighted)
    sprintf('Mean RMSE multilead Kalman filter (%2.2f)\n',rmse_multikf)
    
    
    % Fusioned
    
    
    % Plots
    h = figure;
    ax(1)=subplot(3,1,1);
    plot(fhr(:,chans),'*')
    hold on
    set(gca,'ColorOrderIndex',1)
    plot(fhr_estim)
    plot(fhr_ref,'m','LineWidth',2)
    plot(fhr_final,'--k','LineWidth',2)
    ylim([100 200])
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
win0 = 5;   % window size for SQI calculation [s]
olap0 = 1;  % overalaping [s]
NCHAN = Data.numChan;
win = win0*fs;
olap = olap0*fs;

%% Obtaining segment-wise FHR, references and
k = 1;
seg = 1;
fhr = zeros(max(Data.sqi.segment),max(Data.sqi.channel));

fhr_ref = zeros(max(Data.sqi.segment),1);
fqrs = cell(size(fecg,1),1);
for j = 1:size(fecg,1)
   fqrs{j} = OSET_MaxSearch(fecg(j,:),1.3/fs); % 2.3 Hz is the expected FHR parameter
end


while k < Data.length-win       % Loop through segments
    fref = qrs_ref(qrs_ref>k&qrs_ref<k+win)-k;
    for ch = 1:NCHAN                 % Loop through channels        
        fqrsloc = fqrs{ch}(fqrs{ch}>=k&fqrs{ch}<k+win)-k+1;
        fhr(seg,ch) = fs*60/mean(diff(fqrsloc));  % in bpm
    end
    %sqi(seg,:) = table2array(Data.sqi(idx,11)); % bsqi13
    fhr_ref(seg,1) =  fs*60/mean(diff(fref));
    k = k+olap;  % move window
    seg = seg+1; % loop table indexes
    disp(['Sample: ' num2str(k)])
end
% using sigmoid function to fit SQI "classes" into regression curve

end
