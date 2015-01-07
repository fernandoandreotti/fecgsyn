clear
clc


path_orig = '/media/fernando/FetalEKG/tuning stuff/';
path_res = '/media/fernando/FetalEKG/tuning stuff/extracted3Hz/';
cd(path_res)
fls_res = dir('*.mat'); % looking for .mat (creating index)
fls_res = arrayfun(@(x)x.name,fls_res,'UniformOutput',false);
cd(path_orig)
slashchar = char('/'*isunix + '\'*(~isunix));
fls_orig = dir('*.mat'); % looking for .mat (creating index)
fls_orig = arrayfun(@(x)x.name,fls_orig,'UniformOutput',false);

for i = 6:length(fls_orig)
    %% Resampling everything
    fs_new = 250;
    load([path_orig fls_orig{i}])
    ch = [1 8 11 22 25 32]; % using 6 channels (decided considering Exp. 1)
    fref = round(out.fqrs{1}/(out.param.fs/fs_new));
    mref = round(out.mqrs/(out.param.fs/fs_new));
    if isempty(out.noise)
        noise = zeros(size(out.mecg));
    else
        noise = sum(cat(3,out.noise{:}),3);
    end
    HF_CUT = 100; % high cut frequency
    LF_CUT = 3; % low cut frequency
    [b_lp,a_lp] = butter(5,HF_CUT/(fs_new/2),'low');
    [b_bas,a_bas] = butter(3,LF_CUT/(fs_new/2),'high');
    fs = out.param.fs;
    orig_fecg = sum(cat(3,out.fecg{:}),3);
    orig_mecg = double(out.mecg);
    mixture = orig_mecg + orig_fecg + noise;     % re-creating abdominal mixture
    mixture = mixture./3000;    % removing gain given during int conversion
    for a = 1:length(ch)
        k = ch(a);
        % mix
        ppmixture(a,:) = resample(mixture(k,:),fs_new,fs);    % reducing number of channels
        lpmix = filtfilt(b_lp,a_lp,ppmixture(a,:));
        ppmixture(a,:) = filtfilt(b_bas,a_bas,lpmix);
        % noise
        ppnoise(a,:) = resample(noise(k,:),fs_new,fs);    % reducing number of channels
        lpmix = filtfilt(b_lp,a_lp,ppnoise(a,:));
        ppnoise(a,:) = filtfilt(b_bas,a_bas,lpmix);
        % fecg
        ppfecg(a,:) = resample(orig_fecg(k,:),fs_new,fs);    % reducing number of channels
        lpmix = filtfilt(b_lp,a_lp,ppfecg(a,:));
        ppfecg(a,:) = filtfilt(b_bas,a_bas,lpmix);
        % mecg
        ppmecg(a,:) = resample(orig_mecg(k,:),fs_new,fs);    % reducing number of channels
        lpmix = filtfilt(b_lp,a_lp,ppmecg(a,:));
        ppmecg(a,:) = filtfilt(b_bas,a_bas,lpmix);
    end
    
    % = Reference variables
    % mref
    % fref
    % ppnoise
    % ppfecg
    % ppmecg
    MHR = 60; %     [in bpm]
    FHR = 120; %    [in bpm]
    mbeats = 60*fs*length(mref)/length(ppfecg); % now im bpm
    fbeats = 60*fs*length(fref)/length(ppfecg); % now im bpm
    
    Pmat = sqrt(sum(ppmecg.^2,2)); % maternal power (in each channel)
    Pmat = Pmat.*(MHR/mbeats);                          % normalized power
    Pfet = sqrt(sum(ppfecg.^2,2)); % fetal power (in each channel)
    Pfet = Pfet.*(FHR/fbeats);                          % normalized power
    Pn   = sqrt(sum(ppnoise.^2,2)); % noise power (in each channel)
    
    SNRfm = 10*log10(Pfet./Pmat);
    SNRnf = 10*log10(Pn./Pfet);
    
    files = cellfun(@(x) ~isempty(x), regexp(fls_res,['rec' num2str(i)]));
    files = find(files);
    for x=1:length(files)
        load([path_res fls_res{files(x)}])
        % = Test variables
        % residuals
        [SNRfm,SNRnf] = calcSNR(residual,mref,fref,fs_new);
        
        %% Calculate SNR
        % = Summarized metric
        
    end
    
end
