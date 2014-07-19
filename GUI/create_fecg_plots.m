function [ f_handles , noise_handles, mecg_handle] = create_fecg_plots( out , CH_CANC)
%CREATE_FECG_PLOTS Creates plots based on out_struct for the GUI
%   
disp('Generating plots ..')

if nargin < 3
    CH_CANC = 5;
end

f_handles = [];
noise_handles = [];
param = out.param;

%% Standard plots preparation

col = [1,0,0; % red
    0,0,1; % blue
    0,0.8,0; % green
    0.4,0.4,0; % dark yellow
    0,0.8,0.8; % cyan
    0.4,0,0.8; % dark magenta
    0.8,0.4,1; % light magenta
    0.4,0.4,1]; % lilac

NB_FOETUSES = size(out.f_model,1); % number of foetuses figured out from the number of foetal hea

LINE_WIDTH = 3;
FONT_SIZE = 15;
NB_EL2PLOT = 3; % number of electodes to plot
PACE = 2;



%% Standard plot: Some generated AECG

% == plots a few final AECG channels
tmp_handle = figure('name','Some generated AECG');
set(tmp_handle, 'Visible', 'off');
f_handles = [f_handles, tmp_handle];

    if NB_EL2PLOT<NB_EL2PLOT*PACE
        compt = 0;
        tm = 1/out.param.fs:1/out.param.fs:out.param.n/out.param.fs;
        ax = zeros(NB_EL2PLOT,1);
        for ee=1:PACE:PACE*NB_EL2PLOT
            compt = compt+1;
            ax(compt) = subplot(NB_EL2PLOT,1,compt);plot(tm,out.mixture(ee,:),...
                'color',col(compt,:),'LineWidth',LINE_WIDTH);
            xlabel('Time [sec]'); ylabel('Amplitude [NU]');
            set(gca,'FontSize',FONT_SIZE);
        end
        set(findall(gcf,'type','text'),'fontSize',FONT_SIZE);
        linkaxes(ax,'x');
    else
        error('not enough input channel for plotting with default configuration \n');
    end

%% Standard plot: VCG plots
% == plot mother and foetuse(s) VCGs
tmp_handle = figure('name','VCG plots');
set(tmp_handle, 'Visible', 'off');
f_handles = [f_handles, tmp_handle];

    LegCell = cell(NB_FOETUSES*2,1);
    for vv=1:3
        % == plot
        ax(2*vv-1)=subplot(3,2,2*vv-1); plot(tm,out.m_model.VCG(vv,:),'color',col(2*vv-1,:),'LineWidth',LINE_WIDTH);
        hold on, plot(tm(out.mqrs),out.m_model.VCG(vv,out.mqrs),'+k','LineWidth',LINE_WIDTH);
        legend(['mother VCG channel: ' int2str(vv)],'MQRS');
        set(gca,'FontSize',FONT_SIZE);
        xlabel('Time [sec]'); ylabel('Amplitude [NU]');
        for fet=1:NB_FOETUSES
            ax(2*vv)=subplot(3,2,2*vv); plot(tm,out.f_model{fet}.VCG(vv,:),'color',col(2*vv+fet,:),'LineWidth',LINE_WIDTH);
            hold on, plot(tm(out.fqrs{fet}),out.f_model{fet}.VCG(vv,out.fqrs{fet}),'+k','LineWidth',LINE_WIDTH);
            LegCell(2*(fet-1)+1) = {['foetus ' int2str(fet) ' VCG channel: ' int2str(vv)]};
            LegCell(2*fet) = {['FQRS ' int2str(fet)]};
        end
        legend(LegCell);
        set(gca,'FontSize',FONT_SIZE);
        xlabel('Time [sec]'); ylabel('Amplitude [NU]');
    end
    linkaxes(ax,'x'); xlim([0 tm(end)]);
    set(findall(gcf,'type','text'),'fontSize',FONT_SIZE);


%% Standard plot: Projected FECG and MECG before being mixed
% == plot the projection of mother and foetuse(s) VCGs
tmp_handle = figure('name','Projected FECG and MECG before being mixed');
set(tmp_handle, 'Visible', 'off');
f_handles = [f_handles, tmp_handle];

    LegCell = cell(NB_FOETUSES+1,1);
    GAIN_F = 1;
    if NB_EL2PLOT<NB_EL2PLOT*PACE
        compt = 0;
        tm = 1/out.param.fs:1/out.param.fs:out.param.n/out.param.fs;
        ax = zeros(NB_EL2PLOT,1);
        for ee=1:PACE:PACE*NB_EL2PLOT
            compt = compt+1;
            ax(compt) = subplot(NB_EL2PLOT,1,compt);plot(tm,out.mecg(ee,:),...
                'color','b','LineWidth',LINE_WIDTH);
            LegCell(1) = {'MECG'};
            for fet=1:NB_FOETUSES
                hold on, plot(tm,GAIN_F*out.fecg{fet}(ee,:),'color',col(fet+3,:),'LineWidth',LINE_WIDTH);
                LegCell(1+fet) = {['FECG ' int2str(fet) ', Gain: ' int2str(GAIN_F)]};
            end
            xlabel('Time [sec]'); ylabel('Amplitude [NU]');
            set(gca,'FontSize',FONT_SIZE);
            legend(LegCell);
        end
        set(findall(gcf,'type','text'),'fontSize',FONT_SIZE);
        linkaxes(ax,'x');
    else
        error('not enough input channel for plotting with default configuration \n');
    end


%% Standard plot: Volume conductor
tmp_handle = create_electrode_plot(param, out);

f_handles = [f_handles, tmp_handle];

%% Standard plot: Heart rate
% == plot mother and foetus heart rates
tmp_handle = figure('name','Heart rate');
set(tmp_handle, 'Visible', 'off');
f_handles = [f_handles, tmp_handle];

    fhr = cell(NB_FOETUSES,1);
    legstr = cell(NB_FOETUSES+1,1);
    mhr = 60./diff(out.mqrs/1000);
    plot(tm(out.mqrs(1:end-1)),mhr,'color','r','LineWidth',LINE_WIDTH,'LineStyle','--');
    legstr{1} = 'MQRS';
    for ff=1:NB_FOETUSES
        fhr{ff} = 60./diff(out.fqrs{ff}/1000);
        hold on, plot(tm(out.fqrs{ff}(1:end-1)),fhr{ff},'color',col(ff+1,:),'LineWidth',LINE_WIDTH,'LineStyle','-');
        legstr{ff+1} = ['FQRS' int2str(ff)];
    end
    legend(legstr);
    xlabel('Time [sec]'); ylabel('FHR [bpm]')
    set(gca,'FontSize',FONT_SIZE);
    set(findall(gcf,'type','text'),'fontSize',FONT_SIZE);
    legend boxoff;







%% Noise plots preparation
noise_handles = [];

% % == constants
% fs = param.fs;
% N = param.n;
% AR_ORDER = 12; % number of poles
% FS_NSTDB = 360; % sampling frequency of NSTDB
% LG_NSTDB = FS_NSTDB*29; % number of points in NSTDB
% NP_NSTDB = 20*FS_NSTDB; % number of points to select in NSTDB records to generate the AR coefficients
% N_SAMP = floor(N/(fs/FS_NSTDB)); % N samples at fs correspond to N_SAMP at FS_NSTDB
% NB_EL = size(epos,1); % number of electrodes
% a = out.noise_misc.a;
% 
% 
% %% Noise plot: Poles before and after being shifted
% tmp_handle = figure('name','Poles before and after being shifted');
% set(tmp_handle, 'Visible', 'off');
% noise_handles = [noise_handles, tmp_handle];
% 
% % TODO: add 'a' to the out struct and then use it here
% [~,hp_1,~] = zplane(1,roots([1 a(:,1)'])); 
% set(hp_1,'color','b','LineWidth',3,'MarkerSize',10);
% hold on, [~,hp_2,~] = zplane(1,roots([1 a(:,end)']));
% set(hp_2,'color','r','LineWidth',3,'MarkerSize',10,'Marker','o');
% 
% % set(gca,'FontSize',FONT_SIZE);
% % set(findall(gcf,'type','text'),'fontSize',FONT_SIZE); 
% xlim([-1 1]); ylim([-1 1]);
% 
% 
% 
% 
% %% Noise plot: AR+PCA generated noise
% 
% % == plot the noise generate using AR model and PCA
%     tmp_handle = figure('name','AR+PCA generated noise');
%     set(tmp_handle, 'Visible', 'off');
%     noise_handles = [noise_handles, tmp_handle];
%     
%     tm = 1/fs:1/fs:N/fs;
%     for cc=1:3
%         subplot(3,1,cc); plot(tm,noise_ar(:,cc),'color',col(cc,:),'LineWidth',LINE_WIDTH);
%         xlim([0 10]);
%         set(gca,'FontSize',FONT_SIZE);
%         set(findall(gcf,'type','text'),'fontSize',FONT_SIZE); 
%     end
%     xlabel('Time [sec]'); ylabel('Amplitude [NU]');
%     set(gca,'FontSize',FONT_SIZE);
%     set(findall(gcf,'type','text'),'fontSize',FONT_SIZE);  
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %% Generated ECG plot
% % generate_ecg_mixture()
% 
% 
% 
% FONT_SIZE = 15;
%    LINE_WIDTH = 2;
%    NB_FET = length(fecg);
%    NB_NOISE = length(noise);
%    f_handle = figure;
%    if debug == 11 % corresponds to running the code from the gui
%        set(f_handle, 'Visible', 'off')
%    end
%    % plot maternal signal
%    subplot(NB_FET+NB_NOISE+2,1,1)
%    plot(mecg(1,:),'k','LineWidth',LINE_WIDTH)
%    title('Maternal signal ch1')
%    set(gca,'FontSize',FONT_SIZE)
%    % plot fetal signals
%    for i=1:NB_FET
%        subplot(NB_FET+NB_NOISE+2,1,1+i)
%        plot(fecg{i}(1,:),'b','LineWidth',LINE_WIDTH)
%        title(['Fetal signal ' num2str(i) ' ch1'])
%        set(gca,'FontSize',FONT_SIZE)
%    end
%    % plot noises
%    for i = 1:NB_NOISE
%        subplot(NB_FET+NB_NOISE+2,1,1+NB_FET+i)
%        plot(noise{i}(1,:),'r','LineWidth',LINE_WIDTH)
%        title(['Noise signal ' num2str(i) ' ch1'])
%        set(gca,'FontSize',FONT_SIZE)
%    end
%    % plot resulting signal
%    subplot(NB_FET+NB_NOISE+2,1,NB_FET+NB_NOISE+2)
%    plot(mixture(1,:),'m','LineWidth',LINE_WIDTH)
%    title('Resulting signal for ch1')
%    set(gca,'FontSize',FONT_SIZE);
%    xlabel('Sample Number','FontSize',FONT_SIZE)
%    set(findall(gcf,'type','text'),'fontSize',FONT_SIZE);   
% 



%% MECG Preparation

fs = 1000;
param = out.param;
peaks = adjust_mqrs_location(out.mixture(CH_CANC,:),out.mqrs,param.fs,0);
[residual, ~] = mecg_cancellation(peaks,out.mixture(CH_CANC,:),'TS-CERUTTI',20,2,1000,11);
tm = 1/fs:1/fs:length(residual)/fs;
ecg = out.mixture(CH_CANC,:);
%% MECG Plot: MECG cancellation
mecg_handle = figure('name','MECG cancellation');
set(mecg_handle, 'Visible', 'off');

plot(tm,ecg,'LineWidth',3);
hold on, plot(tm,ecg-residual,'--k','LineWidth',3);
hold on, plot(tm,residual-1.5,'--r','LineWidth',3);
hold on, plot(tm(peaks),ecg(peaks),'+r','LineWidth',2);
hold off
legend('mixture','template','residual','MQRS'); 
title('Template subtraction for extracting the FECG');
xlabel('Time [sec]'); ylabel('Amplitude [NU]')






end

