function FECGSYN_genresults(path_orig,fs,ch,exp3,debug)
% function FECGSYN_genresults(path_orig,fs,ch,exp3,debug)
% 
% this script generates a series of abdominal mixtures, containing i) a
% stationary case and ii) non-stationary case (when adding breathing
% effects, foetal movement etc).
% Experiment 2 - Statistics (F1,MAE,PPV,SE)
% Experiment 3 - Morphologycal Analysis
%
% Input:
% path_orig:        Path for original dataset
% fs:               Sampling frequency
% ch:               Channels to be used
% exp3:            Boolean, if 0 runs exp2 and 1 exp3
%
% Examples:
% TODO
%
% See also:
% TODO
% 
% --
% fecgsyn toolbox, version 1.2, Jan 2017
% Released under the GNU General Public License
%
% Copyright (C) 2014  Joachim Behar & Fernando Andreotti
% University of Oxford, Intelligent Patient Monitoring Group - Oxford 2014
% joachim.behar@oxfordalumni.org, fernando.andreotti@eng.ox.ac.uk
%
% 
% For more information visit: https://www.physionet.org/physiotools/ipmcode/fecgsyn/
% 
% Referencing this work
%
%   Behar Joachim, Andreotti Fernando, Zaunseder Sebastian, Li Qiao, Oster Julien, Clifford Gari D. 
%   An ECG simulator for generating maternal-foetal activity mixtures on abdominal ECG recordings. 
%   Physiological Measurement.35 1537-1550. 2014.
%
% Last updated : 10-03-2016
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

slashchar = char('/'*isunix + '\'*(~isunix));
if debug,  mkdir([path_orig 'plots' slashchar]), end


%% Experiment 1
% stat = stats_struct;
% cc = linspecer(5);
% ch = [2,4,6,8,12,16];
% FONT_SIZE = 14;
% LWIDTH = 2;
% 
% for kk=1:length(stat)
%     mean_FASTICA_DEF(kk) = 100.*mean(stat{kk}.stats_FASTICA_DEF(:,1));
%     median_FASTICA_DEF(kk) = 100.*median(stat{kk}.stats_FASTICA_DEF(:,1));
%     mean_FASTICA_SYM(kk) = 100.*mean(stat{kk}.stats_FASTICA_SYM(:,1));
%     median_FASTICA_SYM(kk) = 100.*median(stat{kk}.stats_FASTICA_SYM(:,1));
%     mean_JADEICA(kk) = 100.*mean(stat{kk}.stats_JADEICA(:,1));
%     median_JADEICA(kk) = 100.*median(stat{kk}.stats_JADEICA(:,1));
%     mean_pca(kk) = 100.*mean(stat{kk}.stats_pca(:,1));
%     median_pca(kk) = 100.*median(stat{kk}.stats_pca(:,1));
% end
% 
% for kk=1:length(stat)
%     FASTICA_DEF(:,kk) = 100.*stat{kk}.stats_FASTICA_DEF(:,1);
%     FASTICA_SYM(:,kk) = 100.*stat{kk}.stats_FASTICA_SYM(:,1);
%     JADEICA(:,kk) = 100.*stat{kk}.stats_JADEICA(:,1);
%     pca(:,kk) = 100.*stat{kk}.stats_pca(:,1);
% end
% 
% figure(1)
% hold on
% plot(ch,mean_JADEICA,'s--','linewidth',LWIDTH,'MarkerFaceColor',cc(1,:),'Color',cc(1,:),'linewidth',LWIDTH)
% plot(ch,median_JADEICA,'s-','linewidth',LWIDTH,'MarkerFaceColor',cc(1,:),'Color',cc(1,:),'linewidth',LWIDTH)
% plot(ch,mean_FASTICA_SYM,'d--','linewidth',LWIDTH,'MarkerFaceColor',cc(2,:),'Color',cc(2,:),'linewidth',LWIDTH)
% plot(ch,median_FASTICA_SYM,'d-','linewidth',LWIDTH,'MarkerFaceColor',cc(2,:),'Color',cc(2,:),'linewidth',LWIDTH)
% plot(ch,mean_FASTICA_DEF,'s--','linewidth',LWIDTH,'MarkerFaceColor',cc(3,:),'Color',cc(3,:),'linewidth',LWIDTH)
% plot(ch,median_FASTICA_DEF,'s-','linewidth',LWIDTH,'MarkerFaceColor',cc(3,:),'Color',cc(3,:),'linewidth',LWIDTH)
% plot(ch,mean_pca,'o--','linewidth',LWIDTH,'MarkerFaceColor',cc(4,:),'Color',cc(4,:),'linewidth',LWIDTH)
% plot(ch,median_pca,'o-','linewidth',LWIDTH,'MarkerFaceColor',cc(4,:),'Color',cc(4,:),'linewidth',LWIDTH)
% xlim([2,32])
% legend('mean ICA (JADE)','median ICA (JADE)','mean ICA (sym. FAST-ICA)','median ICA (sym. FAST-ICA)',...
%     'mean ICA (defl. FAST-ICA)','median ICA (defl. FAST-ICA)','mean PCA','median PCA')
% set(findall(gcf,'type','text'),'fontSize',FONT_SIZE);
% set(gca,'FontSize',FONT_SIZE)
% set(gca,'XTick',ch)
% ylabel('F_1 (in %)')
% xlabel('Number of Channels')
% box on
% % title('PCA reduction')
% hold off


%% Experiment 2
% == Parameters
INTERV = round(0.05*fs); % BxB acceptance interval
TEMP_SAMPS = round(60*fs); % samples used for building templates
fs_new = 250;

% Run through extracted datasets
cd(path_orig)

% abdominal mixtures
fls_orig = dir([path_orig '*.mat']); % looking for .mat (creating index)
fls_orig = arrayfun(@(x)x.name,fls_orig,'UniformOutput',false);
% extracted files
if ~exp3
    path_ext = [path_orig 'exp2' slashchar];
    fls_ext = dir([path_ext '*.mat']); % looking for .mat (creating index)
else
    path_ext = [path_orig 'exp3' slashchar];
    fls_ext = dir([path_ext '*.mat']); % looking for .mat (creating index)
    % Preprocessing Filter coefficiens
    % high-pass filter
    Fstop = 0.5;  % Stopband Frequency
    Fpass = 1;    % Passband Frequency
    Astop = 20;   % Stopband Attenuation (dB)
    Apass = 0.1;  % Passband Ripple (dB)
    h = fdesign.highpass('fst,fp,ast,ap', Fstop, Fpass, Astop, Apass, fs_new);
    Hhp = design(h, 'butter', ...
        'MatchExactly', 'stopband', ...
        'SOSScaleNorm', 'Linf', ...
        'SystemObject', true);
    [b_hp,a_hp] = tf(Hhp);
    % low-pass filter
    Fpass = 90;   % Passband Frequency
    Fstop = 100;  % Stopband Frequency
    Astop = 20;   % Stopband Attenuation (dB)
    Apass = 0.1;    % Passband Ripple (dB)
    h = fdesign.lowpass('fp,fst,ap,ast', Fpass, Fstop, Apass, Astop, fs_new);
    Hlp = design(h, 'butter', ...
        'MatchExactly', 'stopband', ...
        'SOSScaleNorm', 'Linf');
    [b_lp,a_lp] = tf(Hlp);
    clear Fstop Fpass Astop Apass h Hhp Hlp
    
end
fls_ext = arrayfun(@(x)x.name,fls_ext,'UniformOutput',false);
idx = cellfun(@(x) strcmp(x(1:3),'rec'),fls_ext); % ignoring files not begining with 'rec'
fls_ext(~idx) = [];

% pre-allocation
stats.JADEICA = zeros(length(fls_orig),4);
stats.PCA = zeros(length(fls_orig),4);
stats.tsc = zeros(length(fls_orig),4);
stats.tspca = zeros(length(fls_orig),4);
stats.tsekf = zeros(length(fls_orig),4);
stats.alms = zeros(length(fls_orig),4);
stats.arls = zeros(length(fls_orig),4);
stats.aesn = zeros(length(fls_orig),4);
morph.JADEICA = cell(length(fls_orig),7);
morph.PCA = cell(length(fls_orig),7);
morph.tsc = cell(length(fls_orig),7);
morph.tspca = cell(length(fls_orig),7);
morph.tsekf = cell(length(fls_orig),7);
morph.alms = cell(length(fls_orig),7);
morph.arls = cell(length(fls_orig),7);
morph.aesn = cell(length(fls_orig),7);


%% Check specific results

% % for ff =  [50    337    547    582    589]
% %     file = strcat(path_ext,['rec' num2str(ff) '_JADEICA.mat']);
% %     load(file)
% %     file = strcat(path_orig,fls_orig{ff});
% %     load(file)
% %     fecg = double(out.fecg{1}(ch,:)); % selecting channels
% %     %= Resampling original data to match extracted (fs - if necessary)
% %     if size(out.mecg,2) ~= size(residual,2)
% %         % fref:          fetal QRS reference
% %         % fecgref:       fetal ECG signals (before mixture)
% %         fref = floor(out.fqrs{1}.*(size(residual,2)/size(out.mecg,2)));
% %         fecgref = zeros(length(ch),size(residual,2));
% %         for k = 1:size(fecgref,1)
% %             fecgref(k,:) = resample(fecg(k,:),fs,out.param.fs);
% %         end
% %     else
% %         fecgref = fecg;
% %         fref = out.fqrs{1};
% %     end
% %     ax(1)=subplot(1,2,1)
% %     fecgref = (fecgref'*diag(1./max(fecgref')))';
% %     plot(fecgref(4,:),'k')
% %     hold on
% %     for lres = 1:5
% %         res((lres-1)*15000+1:lres*15000) = residual(maxch(lres),(lres-1)*15000+1:lres*15000);
% %     end
% %     plot(res,'LineWidth',1.5)
% %     plot(fref,1,'or')
% %     title('Baseline')
% %     xlabel(['#components = ' num2str(size(residual,1))])
% %     
% %     
% %     file = strcat(path_ext,['rec' num2str(ff+1) '_JADEICA.mat']);
% %     load(file)
% %     file = strcat(path_orig,fls_orig{ff+1});
% %     load(file)
% %     fecg = double(out.fecg{1}(ch,:)); % selecting channels
% %     %= Resampling original data to match extracted (fs - if necessary)
% %     if size(out.mecg,2) ~= size(residual,2)
% %         % fref:          fetal QRS reference
% %         % fecgref:       fetal ECG signals (before mixture)
% %         fref = floor(out.fqrs{1}.*(size(residual,2)/size(out.mecg,2)));
% %         fecgref = zeros(length(ch),size(residual,2));
% %         for k = 1:size(fecgref,1)
% %             fecgref(k,:) = resample(fecg(k,:),fs,out.param.fs);
% %         end
% %     else
% %         fecgref = fecg;
% %         fref = out.fqrs{1};
% %     end
% %     ax(2)=subplot(1,2,2)
% %     fecgref = (fecgref'*diag(1./max(fecgref')))';
% %     plot(fecgref(4,:),'k')
% %     hold on
% %     for lres = 1:5
% %         res((lres-1)*15000+1:lres*15000) = residual(maxch(lres),(lres-1)*15000+1:lres*15000);
% %     end
% %     plot(res,'LineWidth',1.5)
% %     plot(fref,1,'or')
% %     linkaxes(ax,'x')        
% %     title('Case 0')
% %     xlabel(['#components = ' num2str(size(residual,1))])
% %     
% %     savefig(['plot_' num2str(ff)])
% %     close
% % end



% = Runs through list of extracted files
for i = filesproc%length(fls_ext)
    % for i = randperm(length(fls_ext))
    disp(fls_ext{i})
    fprintf('Data %d out of %d \n',i,length(fls_ext));
    
    %= loading extracted file
    [rec,met] = strtok(fls_ext(i),'_');
    % Figuring out which extraction method was used, possibilities are:
    % (JADEICA,PCA,tsc,tspca,tsekf,alms,arls,aesn)
    method = met{:}(2:end-4);
    if ~strcmp(method,'JADEICA')
        continue
    end
    file = strcat(path_ext,fls_ext(i));
    load(file{:})
    %= loading original file
    origrec = str2double(rec{:}(4:end));
    file = strcat(path_orig,fls_orig(origrec));
    cas = regexp(file{:},'_c[0-7]','match'); % find out which case it depicts
    if isempty(cas)
        cas = {'bas'};
    end
    
    load(file{:});
    fecg = double(out.fecg{1}(ch,:)); % selecting channels
    %= Resampling original data to match extracted (fs - if necessary)
    if size(out.mecg,2) ~= size(residual,2)
        % fref:          fetal QRS reference
        % fecgref:       fetal ECG signals (before mixture)
        fref = floor(out.fqrs{1}.*(size(residual,2)/size(out.mecg,2)));
        fecgref = zeros(length(ch),size(residual,2));
        for k = 1:size(fecgref,1)
            fecgref(k,:) = resample(fecg(k,:),fs,out.param.fs);
        end
    else
        fecgref = fecg;
        fref = out.fqrs{1};
    end
    [elif,~]=strtok(file{:}(end:-1:1),slashchar);
    disp(elif(end:-1:1))
    clear fecg outdata rec file elif k
    
    %= Getting statistics (exp 2)
    if ~exp3
        [F1,MAE,PPV,SE] = Bxb_compare(fref,fqrs,INTERV);
        MAE = MAE*1000/fs_new;
        stats.(method)(origrec,:) = [F1,MAE,PPV,SE]; % dynamic naming
    else %= Getting statistics (exp 3)
        
    end
    clear fecg residual fqrs F1 MAE PPV SE qt_err theight_err outputs
    cd ..
end

save([path_orig 'wksp' num2str(filesproc(end))])

    
    
    
end

end
