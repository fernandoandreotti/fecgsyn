function run_KFHR

% Params
R0 = 1e-3;
Q0 = 1;
p = 1;
win0 = 5;   % window size for SQI calculation [s]
olap0 = 1;  % overalaping [s]
tol = 10;  % hdr tolerance [bpm]
thres = 0; % threshold for SQI to be considered [0-1]

% Loop through
fls = dir('*.mat');
fls = arrayfun(@(x) x.name,fls,'UniformOutput',false);
load('nbmodel.mat')          % decision tree classificator

rmse_final = zeros(length(fls));
rmse_pre1 = rmse_final; rmse_pre2 = rmse_final; rmse_pre3 = rmse_final;
rmse1 = rmse_final; rmse2 = rmse_final; rmse3 = rmse_final;


hdr_final = rmse_final;
hdr_pre1 = hdr_final; hdr_pre2 = hdr_final; hdr_pre3 = hdr_final;
hdr1 = hdr_final; hdr2 = hdr_final; hdr3 = hdr_final;


chans{1} = 5:7;
chans{2} = 1:4;
chans{3} = 1:7;
chans{4} = [1,7];
tabhdr = cell(2,1);
tabhdr{1} = zeros(length(fls),6);
tabhdr{2} = zeros(length(fls),6);
tabhdr{3} = zeros(length(fls),6);
tabrmse = tabhdr;

for c = 4
    % c=3;
    count = 1;
    hdr = [];
    % for thres = 0:0.1:0.9
    for i = 1:length(fls)
        load(fls{i})
        disp(num2str(i))
        %         [sqinorm,fhr,fhr_ref] = sqi_extracted_consensus(Data); % get sqis and references
        %% Run Kalman filter
        %     fhr_ref - FHR reference
        %     fhr     - estimated FHR for each channel
        %     sqi     - SQI for each channel
        [fhr_final,fhr_estim,innov_estim] = multichan_hr_KF(fhr(:,chans{c}),sqinorm(:,chans{c}),p,Q0,R0); % multi and single channel results
        fhr_weight=sum(fhr(:,chans{c}).*sqinorm(:,chans{c})./repmat(sum(sqinorm(:,chans{c}),2),1,length(chans{c})),2); % weighted average on SQI values
        
        % removing bad quality segments!!
        rem = isnan(fhr_ref)';
        if any(rem)
            fhr_ref(rem) = []; fhr_estim(rem,:) = [];  fhr_final(rem,:) = []; fhr(rem,:) = []; fhr_weight(rem) = [];
        end
        
        % HDR results
        hdr_pre =  1-sum(abs(repmat(fhr_ref,1,length(chans{c}))-fhr(:,chans{c}))>tol)./length(fhr_ref);
        hdr_kfsingle =  1-sum(abs(repmat(fhr_ref,1,length(chans{c}))-fhr_estim)>tol)./length(fhr_ref);
        hdr_weighted = 1-sum(abs(fhr_ref-fhr_weight)>tol)/length(fhr_ref);
        hdr_multikf = 1-sum(abs(fhr_ref-fhr_final)>tol)/length(fhr_ref);
        tabhdr{c}(i,1) = mean(hdr_pre);
        tabhdr{c}(i,2) = mean(hdr_kfsingle);
        tabhdr{c}(i,3) = hdr_weighted;
        tabhdr{c}(i,4) = hdr_multikf;
        tabhdr{c}(i,5) = max(hdr_pre);
        tabhdr{c}(i,6) = max(hdr_kfsingle);
        % RMSE results
        
        rmse_pre =  sqrt(sum((repmat(fhr_ref,1,length(chans{c}))-fhr(:,chans{c})).^2)/length(fhr_ref));
        rmse_kfsingle =  sqrt(sum((repmat(fhr_ref,1,length(chans{c}))-fhr_estim).^2)/length(fhr_ref));
        rmse_weighted = sqrt(sum((fhr_ref-fhr_weight).^2)/length(fhr_ref));
        rmse_multikf = sqrt(sum((fhr_ref-fhr_final).^2)/length(fhr_ref));
        tabrmse{c}(i,1) = mean(rmse_pre);
        tabrmse{c}(i,2) = mean(rmse_kfsingle);
        tabrmse{c}(i,3) = rmse_weighted;
        tabrmse{c}(i,4) = rmse_multikf;
        tabrmse{c}(i,5) = min(rmse_pre);
        tabrmse{c}(i,6) = min(rmse_kfsingle);
        
        % Fusioned
        
        
        % Plots
                h = figure;
                ax(1)=subplot(3,1,1);
                plot(fhr(:,chans{c}),'*')
                hold on
                set(gca,'ColorOrderIndex',1)
                plot(fhr_estim)
                plot(fhr_ref,'m','LineWidth',2)
                plot(fhr_final,'--k','LineWidth',2)
                ylim([100 200])
                ylabel('FHR (in bpm)')
                ax(2)=subplot(3,1,2);
                set(gca,'ColorOrderIndex',1)
                plot(100.*sqinorm(:,chans{c}))
                ylim([0 100])
                ylabel('SQI (in %)')
                ax(3)=subplot(3,1,3);
                set(gca,'ColorOrderIndex',1)
                plot(innov_estim)
                ylabel('Innovation (\nu)')
                linkaxes(ax,'x')
%                 xlim([990 1050])
        %         saveas(h,sprintf('FHR%d_chan%d.png',i,c)); % will create FIG1, FIG2,...
        close
        %         %
        clear Data fhr fhr_ref fhr_final fhr_estim fhr_weight
    end
    % end
    tab(count,1)=mean(badqual);
    tab(count,2)=mean(tabhdr{3,1}(:,4));
    count = count+1;
    clear badqual
end

end

function [sqinorm,fhr,fhr_ref] = prepare4KFHR(Data)
load('nbmodel.mat')          % decision tree classificator
% Params
win0 = 5;   % window size for SQI calculation [s]
olap0 = 1;  % overalaping [s]
acceptint = 0.050; % 50 ms
coding = [0, 0.4764, 0.8047, 0.9660, 1.0000]';
NCHAN = Data.numChan;
win = win0*Data.sampRate;
olap = olap0*Data.sampRate;
mref_all = Data.mRef.samplestamp;       % maternal QRS reference
fref_all = Data.fRef.samplestamp;       % maternal QRS reference


%% Get SQI decision predictions
features = table2array(Data.sqi(:,5:end));
features(isnan(features)) = 0; % substitute NaNs for zeros
sqinorm = zeros(max(Data.sqi.segment),max(Data.sqi.channel));
for ch = 1:NCHAN                 % Loop through channels
    idx = Data.sqi.channel==ch;
    [scores] = nbmodel.posterior(features(idx,:));   
    sqinorm(:,ch)=scores*coding;
    feats(:,:,ch) = features(idx,:);
end

%% Obtaining segment-wise FHR, references and
k = 1;
seg = 1;
fhr = zeros(max(Data.sqi.segment),max(Data.sqi.channel));
fqrsacc =  fhr;
fhr_ref = zeros(max(Data.sqi.segment),1);
fqrs = QRS_detect(Data,true,'MaxSearch');
while k < Data.length-win       % Loop through segments
    mref = mref_all(mref_all>k&mref_all<k+win)-k;
    fref = fref_all(fref_all>k&fref_all<k+win)-k;
    for ch = 1:NCHAN                 % Loop through channels
        fqrsloc = fqrs{ch}(fqrs{ch}>=k&fqrs{ch}<k+win)-k+1;
        fhr(seg,ch) = Data.sampRate*60/mean(diff(fqrsloc));  % in bpm
        fqrsacc(seg,ch) = bsqi(fref,fqrsloc,acceptint,Data.sampRate);       
    end
    %sqi(seg,:) = table2array(Data.sqi(idx,11)); % bsqi13
    fhr_ref(seg,1) =  Data.sampRate*60/mean(diff(fref));
    k = k+olap;  % move window
    seg = seg+1; % loop table indexes
    disp(['Sample: ' num2str(k)])
end
% using sigmoid function to fit SQI "classes" into regression curve

%% Saving results
save([Data.name '_kfready.mat'],'feats','fqrsacc','fhr','fhr_ref','sqinorm')
end
