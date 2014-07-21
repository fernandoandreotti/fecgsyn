function [ f_handles , noise_handles, gen_ecg_handle, mecg_handle] = FECGSYN_UI_create_fecg_plots( out , choice, CH_CANC)
%CREATE_FECG_PLOTS Creates plots based on out_struct for the GUI
%   choice is a 4-element binary vector corresponding to choosing which of
%   the [standard, noise, gen ecg, mecg] plots to generate. 
%
% Mohsan Alvi (mohsan.alvi@eng.ox.ac.uk) - July 2014
disp('Generating plots ..')

if nargin < 2
    choice = [1,1,1,1];
end
if nargin < 3
    CH_CANC = 5;
end

f_handles = [];
noise_handles = [];
gen_ecg_handle = [];
mecg_handle = [];
param = out.param;
debug = 11;

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

LINE_WIDTH = 2;
FONT_SIZE = 15;
FONT_SIZE_SMALL = 10;
NB_EL2PLOT = 3; % number of electodes to plot
PACE = 2;



%% Standard plot: Abdominal ECG mixture
if choice(1)
    % == plots a few final AECG channels
    tmp_handle = figure('name','Abdominal ECG mixture');
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
                %xlabel('Time [sec]'); ylabel('Amplitude [NU]');
                set(gca,'FontSize',FONT_SIZE);
            end
            xlabel(ax(end), 'Time [sec]'); ylabel(ax(2), 'Amplitude [NU]');
            set(findall(gcf,'type','text'),'fontSize',FONT_SIZE);
            linkaxes(ax,'x');
        else
            error('not enough input channel for plotting with default configuration \n');
        end

    %% Standard plot: Vectorcardiogram
    % == plot mother and foetuse(s) VCGs
    tmp_handle = figure('name','Vectorcardiogram');
    set(tmp_handle, 'Visible', 'off');
    f_handles = [f_handles, tmp_handle];

        LegCell = cell(NB_FOETUSES*2,1);
        for vv=1:3
            % == plot
            ax(2*vv-1)=subplot(3,2,2*vv-1); plot(tm,out.m_model.VCG(vv,:),'color',col(2*vv-1,:),'LineWidth',LINE_WIDTH);
            hold on, plot(tm(out.mqrs),out.m_model.VCG(vv,out.mqrs),'+k','LineWidth',LINE_WIDTH);
            lgnd = legend(['mother VCG channel: ' int2str(vv)],'MQRS');
            set(gca,'FontSize',FONT_SIZE_SMALL);
            %set(lgnd,'FontSize', FONT_SIZE_SMALL);
            %xlabel('Time [sec]'); ylabel('Amplitude [NU]');
            for fet=1:NB_FOETUSES
                ax(2*vv)=subplot(3,2,2*vv); plot(tm,out.f_model{fet}.VCG(vv,:),'color',col(2*vv+fet,:),'LineWidth',LINE_WIDTH);
                hold on, plot(tm(out.fqrs{fet}),out.f_model{fet}.VCG(vv,out.fqrs{fet}),'+k','LineWidth',LINE_WIDTH);
                LegCell(2*(fet-1)+1) = {['foetus ' int2str(fet) ' VCG channel: ' int2str(vv)]};
                LegCell(2*fet) = {['FQRS ' int2str(fet)]};
            end
            legend(LegCell);
            set(gca,'FontSize',FONT_SIZE_SMALL);
            xlabel('Time [sec]'); ylabel('Amplitude [NU]');
        end
        linkaxes(ax,'x'); xlim([0 tm(end)]);
        set(findall(gcf,'type','text'),'fontSize',FONT_SIZE_SMALL);


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
    tmp_handle = FECGSYN_UI_create_electrode_plot(param, out);

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





end

%% Noise plots preparation
% TODO: Check how this will work with more than one noise source
if choice(2)
    % Only generate noise plots if there are noise sources
    if any(strcmp('ntype',fieldnames(param))) && ~isempty(param.ntype)

        % == constants
        a = out.noise_misc.a;
        epos = out.noise_misc.epos;
        noise_ar = out.noise_misc.n_model{1}.VCG';
        ainit = out.noise_misc.ainit;
        
        fs = param.fs;
        N = param.n;
        AR_ORDER = 12; % number of poles
        FS_NSTDB = 360; % sampling frequency of NSTDB
        LG_NSTDB = FS_NSTDB*29; % number of points in NSTDB
        NP_NSTDB = 20*FS_NSTDB; % number of points to select in NSTDB records to generate the AR coefficients
        N_SAMP = floor(N/(fs/FS_NSTDB)); % N samples at fs correspond to N_SAMP at FS_NSTDB
        NB_EL = size(epos,1); % number of electrodes
        


        %% Noise plot: Poles before and after being shifted
        tmp_handle = figure('name','Poles before and after being shifted');
        set(tmp_handle, 'Visible', 'off');
        noise_handles = [noise_handles, tmp_handle];

        % TODO: add 'a' to the out struct and then use it here
        [~,hp_1,~] = zplane(1,roots([1 a(:,1)'])); 
        set(hp_1,'color','b','LineWidth',3,'MarkerSize',10);
        hold on, [~,hp_2,~] = zplane(1,roots([1 a(:,end)']));
        set(hp_2,'color','r','LineWidth',3,'MarkerSize',10,'Marker','o');

        % set(gca,'FontSize',FONT_SIZE);
        % set(findall(gcf,'type','text'),'fontSize',FONT_SIZE); 
        xlim([-1 1]); ylim([-1 1]);




        %% Noise plot: AR+PCA generated noise

        % == plot the noise generate using AR model and PCA
        tmp_handle = figure('name','AR+PCA generated noise');
        set(tmp_handle, 'Visible', 'off');
        noise_handles = [noise_handles, tmp_handle];

        tm = 1/fs:1/fs:N/fs;
        ax = -1*ones(3,1);
        for cc=1:3
            ax(cc) = subplot(3,1,cc); plot(tm,noise_ar(:,cc),'color',col(cc,:),'LineWidth',LINE_WIDTH);
            xlim([0 10]);
            set(gca,'FontSize',FONT_SIZE);
            set(findall(gcf,'type','text'),'fontSize',FONT_SIZE); 
        end
        xlabel(ax(end), 'Time [sec]'); ylabel(ax(2), 'Amplitude [NU]');
        set(gca,'FontSize',FONT_SIZE);
        set(findall(gcf,'type','text'),'fontSize',FONT_SIZE);  




        %% Noise plot: Power Spectral Density plot
        
        % = using the AR coeff computed here
        tmp_handle = figure('name','Power Spectral Density plot');
        set(tmp_handle, 'Visible', 'off');
        noise_handles = [noise_handles, tmp_handle];
        
        [h1,f1] = freqz(1,ainit,512,FS_NSTDB);
        P1 = abs(h1).^2; % power
        P1dB = 10*log10(P1/(mean(P1))); % power in decibels
        plot(f1,P1dB,'LineWidth',LINE_WIDTH);

        [h2,f2] = freqz(1,[1; mean(a,2)],512,FS_NSTDB);
        P2 = abs(h2).^2;
        P2dB = 10*log10(P2/(mean(P2)));
        hold on, plot(f2,P2dB,'--r','LineWidth',LINE_WIDTH);

        xlabel('Frequency [Hz]');
        ylabel('Power [db]');
        set(gca,'FontSize',FONT_SIZE);
        set(findall(gcf,'type','text'),'fontSize',FONT_SIZE);   
        legend('initial AR coefficients','average AR coefficients');
        box off;
        legend boxoff;
        
    end
end


%% Generated ECG plot preparation
% generate_ecg_mixture()
if choice(3)    

    f_model = out.f_model;
    n_model = out.noise_misc.n_model;
    mqrs = out.mqrs;
    fqrs = out.fqrs;
    m_model = out.m_model;
    
    %% Generated ECG plot: 

    [~,~,~,~, gen_ecg_handle] = generate_ecg_mixture(debug,param.SNRfm,...
            param.SNRmn,mqrs,fqrs,param.fs,m_model,f_model{:},n_model{:});
        
end

%% MECG Preparation
if choice(4)
    fs = 1000;
    param = out.param;
    peaks = adjust_mqrs_location(out.mixture(CH_CANC,:),out.mqrs,param.fs,0);
    [residual, ~] = mecg_cancellation(peaks,out.mixture(CH_CANC,:),'TS-CERUTTI',20,2,1000,11);
    tm = 1/fs:1/fs:length(residual)/fs;
    ecg = out.mixture(CH_CANC,:);
    %% MECG Plot: Maternal ECG cancellation
    mecg_handle = figure('name','Maternal ECG cancellation');
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




end

