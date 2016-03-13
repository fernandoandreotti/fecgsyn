function [ f_handles , noise_handles, gen_ecg_handle, mecg_handle] = FECGSYN_UI_create_fecg_plots( out , choice, CH_CANC)
% function [ f_handles , noise_handles, gen_ecg_handle, mecg_handle] = FECGSYN_UI_create_fecg_plots( out , choice, CH_CANC)
% FECGSYN is fruit of the collaboration between the Department of Engineering 
% Science, University of Oxford (DES-OX) and the Institute of Biomedical Engineering, 
% TU Dresden (IBMT-TUD). The authors are Joachim Behar (DES-OX), Fernando Andreotti 
% (IBMT-TUD), Julien Oster (DES-OX), Sebastian Zaunseder (IBMT-TUD) and 
% Gari Clifford (DES-OX). 
%
% The present user interface was contributed by Mohsan Alvi (DES-OX) under
% the supervision of Joachim Behar (DES-OX) and Fernando Andreotti (IBMT-TUD).
%
% --
% fecgsyn toolbox, version 1.1, March 2016
% Released under the GNU General Public License
%
% Copyright (C) 2014  Joachim Behar & Fernando Andreotti
% Oxford university, Intelligent Patient Monitoring Group - Oxford 2014
% joachim.behar@eng.ox.ac.uk, fernando.andreotti@mailbox.tu-dresden.de
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
    CH_CANC = 10;
end

f_handles = {};
noise_handles = {};
gen_ecg_handle = {};
mecg_handle = {};
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

PACE = 1;
%NB_EL = 34; % There are 34 electrodes
NB_EL = size(out.param.elpos, 1); % There are 34 electrodes
NB_EL2PLOT = min([3, NB_EL]); % number of electodes to plot

tm_idx = out.param.n; % plot elements 1:15000 for first 15 seconds
%tm_idx = 60000; % plot elements 1:15000 for first 15 seconds
mqrs_1t15 = out.mqrs(out.mqrs<tm_idx);  % mqrs up to 15 seconds

for i = 1:numel(out.fqrs)
    fqrs_1t15{i} = out.fqrs{i}(out.fqrs{i}<tm_idx);  % fqrs up to 15 seconds for each fetus
end


%% Standard plot: Abdominal ECG mixture
if choice(1)
    % == plots AECG channels
    NB_PLOTS = ceil(NB_EL/3); % 3 electrodes shown on each plot
    
    tmp_handle = struct;
    tmp_handle.title = 'Abdominal ECG mixture';
    tmp_handle.plots = -1 * ones(NB_PLOTS, 1);
    
        tm = 1/out.param.fs:1/out.param.fs:(out.param.n/out.param.fs);

        for pp=1:NB_PLOTS
            if pp == 1 % which electrodes are being plotted
                el = linspace(1, 3, 3);
            else
                el = el + 3;
                el = el(el <= NB_EL);
            end

            tmp_handle.plots(pp) = figure('name', sprintf('Electrodes %d - %d', [el(1), el(end)]));
            set(tmp_handle.plots(pp), 'Visible', 'off');

            compt = 0;
            for ee = 1: length(el) % looping through length of ee to stop correctly at 34th electrode 
                compt = compt+1;
                subplot(3,1,compt);plot(tm(1:tm_idx),out.mixture(el(ee),1:tm_idx),...
                'color',col(compt,:),'LineWidth',LINE_WIDTH);
                xlim([1 4]);
                if ee == 2
                    ylabel('Amplitude [NU]');
                elseif length(el)==1 && ee == 1
                    ylabel('Amplitude [NU]');
                end
                
                if ee ~= length(el)
                    set(gca,'xtick',[])
                    set(gca,'xticklabel',[])
                end
                
                if ee == length(el)
                    xlabel('Time [sec]')
                end
            end
            set(gca,'FontSize',FONT_SIZE_SMALL);
        end
        set(findall(gcf,'type','text'),'fontSize',FONT_SIZE_SMALL);
        %linkaxes(ax,'x');

    

%         if NB_EL2PLOT<NB_EL2PLOT*PACE
%             compt = 0;
%             tm = 1/out.param.fs:1/out.param.fs:out.param.n/out.param.fs;
%             ax = zeros(NB_EL2PLOT,1);
%             for ee=1:PACE:PACE*NB_EL2PLOT
%                 compt = compt+1;
%                 ax(compt) = subplot(NB_EL2PLOT,1,compt);plot(tm,out.mixture(ee,:),...
%                     'color',col(compt,:),'LineWidth',LINE_WIDTH);
%                 %xlabel('Time [sec]'); ylabel('Amplitude [NU]');
%                 set(gca,'FontSize',FONT_SIZE_SMALL);
%             end
%             xlabel(ax(end), 'Time [sec]'); ylabel(ax(2), 'Amplitude [NU]');
%             set(findall(gcf,'type','text'),'fontSize',FONT_SIZE_SMALL);
%             linkaxes(ax,'x');
%         else
%             error('not enough input channel for plotting with default configuration \n');
%         end    
    
    f_handles{end+1} = tmp_handle;
    
    
    %% Standard plot: Vectorcardiogram
    % == plot mother and foetuse(s) VCGs
    tmp_handle = struct;
    tmp_handle.title = 'Vectorcardiogram';
    tmp_handle.plots = -1 * ones(6,1);
%     set(tmp_handle.plots, 'Visible', 'off');

        % Initialising the figures
        for vv = 1:3
            tmp_handle.plots(2*vv-1) = figure('name',['Mother VCG channel ' int2str(vv)]);
            set(tmp_handle.plots(2*vv-1), 'visible', 'off')
            tmp_handle.plots(2*vv) = figure('name', ['Foetus VCG channel: ' int2str(vv)] );
            set(tmp_handle.plots(2*vv), 'visible', 'off')
        end

        LegCell = cell(NB_FOETUSES*2,1);

        for vv=1:3
            % == plot
            figure(tmp_handle.plots(2*vv-1))
            set(tmp_handle.plots(2*vv-1), 'visible', 'off')
            plot(tm(1:tm_idx),out.m_model.VCG(vv,1:tm_idx),'color',col(2*vv-1,:),'LineWidth',LINE_WIDTH);
            hold on, plot(tm(mqrs_1t15),out.m_model.VCG(vv,mqrs_1t15),'+k','LineWidth',LINE_WIDTH);
            lgnd = legend(['mother VCG channel: ' int2str(vv)],'MQRS');
            set(gca,'FontSize',FONT_SIZE_SMALL);
            xlabel('Time [sec]'); ylabel('Amplitude [NU]');
            xlim([1 4]);
            %set(lgnd,'FontSize', FONT_SIZE_SMALL);

            for fet=1:NB_FOETUSES
                figure(tmp_handle.plots(2*vv))
                set(tmp_handle.plots(2*vv), 'visible', 'off')
                plot(tm(1:tm_idx),out.f_model{fet}.VCG(vv,1:tm_idx),'color',col(2*vv+fet,:),'LineWidth',LINE_WIDTH);
                hold on, plot(tm(fqrs_1t15{fet}),out.f_model{fet}.VCG(vv,fqrs_1t15{fet}),'+k','LineWidth',LINE_WIDTH);
                LegCell(2*(fet-1)+1) = {['foetus ' int2str(fet) ' VCG channel: ' int2str(vv)]};
                LegCell(2*fet) = {['FQRS ' int2str(fet)]};
                ylabel('Amplitude [NU]'); xlabel('Time [sec]');
                xlim([1 4]);
            end
            legend(LegCell);
            set(gca,'FontSize',FONT_SIZE_SMALL);
        end
    %         linkaxes(ax,'x'); xlim([0 tm(end)]);
    %         set(findall(gcf,'type','text'),'fontSize',FONT_SIZE_SMALL);

    f_handles{end+1} = tmp_handle;  
    
    
    %% Standard plot: Projected FECG and MECG before being mixed
    % == plot the projection of mother and foetuse(s) VCGs
    NB_PLOTS = ceil(NB_EL/3); % 3 electrodes shown on each plot
    
    tmp_handle = struct;
    tmp_handle.title = 'Projected FECG and MECG before being mixed';
    tmp_handle.plots = -1 * ones(NB_PLOTS, 1);

        for curr_plot = 1:NB_PLOTS
            if curr_plot == 1 % which electrodes are being plotted
                el = linspace(1, 3, 3);
            else
                el = el + 3;
                el = el(el <= NB_EL);
            end
            
            tmp_handle.plots(curr_plot) = figure('name', sprintf('Electrodes %d - %d', [el(1), el(end)]));
            set(tmp_handle.plots(curr_plot), 'Visible', 'off');
            
            LegCell = cell(NB_FOETUSES+1,1);
            GAIN_F = 1;

            compt = 0;

            ax = zeros(NB_EL2PLOT,1);
            for ee= 1: length(el) % looping through length of ee to stop correctly at 34th electrode 
                compt = compt+1;
                ax(compt) = subplot(NB_EL2PLOT,1,compt);plot(tm(1:tm_idx),out.mecg(el(ee),1:tm_idx),...
                    'color','b','LineWidth',LINE_WIDTH);
                LegCell(1) = {'MECG'};
                for fet=1:NB_FOETUSES
                    hold on, plot(tm(1:tm_idx),GAIN_F*out.fecg{fet}(el(ee),1:tm_idx),'color',col(fet+3,:),'LineWidth',LINE_WIDTH);
                    LegCell(1+fet) = {['FECG ' int2str(fet) ', Gain: ' int2str(GAIN_F)]};
                end
                %xlabel('Time [sec]'); ylabel('Amplitude [NU]');
                if ee == 2
                    ylabel('Amplitude [NU]');
                elseif length(el)==1 && ee == 1
                    ylabel('Amplitude [NU]');
                end
                
                if ee ~= length(el)
                    set(gca,'xtick',[])
                    set(gca,'xticklabel',[])
                end
                
                if ee == length(el)
                    xlabel('Time [sec]')
                end
                xlim([1 4]);
                set(gca,'FontSize',FONT_SIZE_SMALL);
                legend(LegCell);
            end
            %set(findall(gcf,'type','text'),'fontSize',FONT_SIZE_SMALL);
            %linkaxes(ax,'x');
        end
%         LegCell = cell(NB_FOETUSES+1,1);
%         GAIN_F = 1;
%         if NB_EL2PLOT<NB_EL2PLOT*PACE
%             compt = 0;
%             tm = 1/out.param.fs:1/out.param.fs:out.param.n/out.param.fs;
%             ax = zeros(NB_EL2PLOT,1);
%             for pp=1:PACE:PACE*NB_EL2PLOT
%                 compt = compt+1;
%                 ax(compt) = subplot(NB_EL2PLOT,1,compt);plot(tm,out.mecg(pp,:),...
%                     'color','b','LineWidth',LINE_WIDTH);
%                 LegCell(1) = {'MECG'};
%                 for fet=1:NB_FOETUSES
%                     hold on, plot(tm,GAIN_F*out.fecg{fet}(pp,:),'color',col(fet+3,:),'LineWidth',LINE_WIDTH);
%                     LegCell(1+fet) = {['FECG ' int2str(fet) ', Gain: ' int2str(GAIN_F)]};
%                 end
%                 xlabel('Time [sec]'); ylabel('Amplitude [NU]');
%                 set(gca,'FontSize',FONT_SIZE_SMALL);
%                 legend(LegCell);
%             end
%             set(findall(gcf,'type','text'),'fontSize',FONT_SIZE_SMALL);
%             linkaxes(ax,'x');
%         else
%             error('not enough input channel for plotting with default configuration \n');
%         end

    f_handles{end+1} = tmp_handle;
    
    
    %% Standard plot: Volume conductor
%     tmp_handle = struct;
%     tmp_handle.title = 'Volume conductor';
%     tmp_handle.plots = FECGSYN_UI_create_electrode_plot(param, out);
%     set(tmp_handle.plots, 'Visible', 'off');
%     f_handles{end+1} = tmp_handle;
    

    %% Standard plot: Heart rate
    % == plot mother and foetus heart rates
    tmp_handle = struct;
    tmp_handle.title = 'Heart rate';
    tmp_handle.plots = figure('name','Heart rate');
    set(tmp_handle.plots, 'Visible', 'off');
    f_handles{end+1} = tmp_handle;

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
        xlabel('Time [sec]'); ylabel('FHR [bpm]');
        xlim([1 4]);
        set(gca,'FontSize',FONT_SIZE_SMALL);
        set(findall(gcf,'type','text'),'fontSize',FONT_SIZE_SMALL);
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
%         NB_EL = size(epos,1); % number of electrodes
        


        %% Noise plot: Poles before and after being shifted
        tmp_handle = struct;
        tmp_handle.title = 'Poles before and after being shifted';
        tmp_handle.plots = figure('name','Poles before and after being shifted');
        set(tmp_handle.plots, 'Visible', 'off');
        noise_handles{end+1} = tmp_handle;

        % TODO: add 'a' to the out struct and then use it here
        [~,hp_1,~] = zplane(1,roots([1 a(:,1)'])); 
        set(hp_1,'color','b','LineWidth',3,'MarkerSize',10);
        hold on, [~,hp_2,~] = zplane(1,roots([1 a(:,end)']));
        set(hp_2,'color','r','LineWidth',3,'MarkerSize',10,'Marker','o');

        % set(gca,'FontSize',FONT_SIZE);
        % set(findall(gcf,'type','text'),'fontSize',FONT_SIZE); 
        xlim([-1 1]); ylim([-1 1]);
    
        set(findall(gcf,'type','text'),'fontSize',FONT_SIZE_SMALL);



        %% Noise plot: AR+PCA generated noise

        % == plot the noise generate using AR model and PCA
        tmp_handle = struct;
        tmp_handle.title = 'AR+PCA generated noise';
        tmp_handle.plots = figure('name','AR+PCA generated noise');
        set(tmp_handle.plots, 'Visible', 'off');
        noise_handles{end+1} = tmp_handle;

        tm = 1/fs:1/fs:N/fs;    % time up to 60s
        
        ax = -1*ones(3,1);
        for cc=1:3
            ax(cc) = subplot(3,1,cc); plot(tm(1:tm_idx),noise_ar(1:tm_idx,cc),'color',col(cc,:),'LineWidth',LINE_WIDTH);
            xlim([1 4]);
            set(gca,'FontSize',FONT_SIZE_SMALL);
            set(findall(gcf,'type','text'),'fontSize',FONT_SIZE_SMALL); 
        end
        xlabel(ax(end), 'Time [sec]'); ylabel(ax(2), 'Amplitude [NU]');
        set(gca,'FontSize',FONT_SIZE_SMALL);
        set(findall(gcf,'type','text'),'fontSize',FONT_SIZE_SMALL);  




        %% Noise plot: Power Spectral Density plot
        
        % = using the AR coeff computed here
        tmp_handle = struct;
        tmp_handle.title = 'Power Spectral Density plot';
        tmp_handle.plots = figure('name','Power Spectral Density plot');
        set(tmp_handle.plots, 'Visible', 'off');
        noise_handles{end+1} = tmp_handle;
        
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
        set(gca,'FontSize',FONT_SIZE_SMALL);
        set(findall(gcf,'type','text'),'fontSize',FONT_SIZE_SMALL);   
        legend('initial AR coefficients','average AR coefficients');
        box off;
        legend boxoff;
        
    end
end


%% Generated ECG plot preparation
% % generate_ecg_mixture()
% if choice(3)    
% 
%     f_model = out.f_model;
%     n_model = out.noise_misc.n_model;
%     mqrs = out.mqrs;
%     fqrs = out.fqrs;
%     m_model = out.m_model;
%     
%     %% Generated ECG plot: 
%         
%     tmp_handle = struct;
%     tmp_handle.title = 'Generated ECG Mixture';
%     [~,~,~,~, tmp_handle.plots] = generate_ecg_mixture(debug,param.SNRfm,...
%             param.SNRmn,mqrs,fqrs,param.fs,m_model,f_model{:},n_model{:});
%     %set(tmp_handle.plots, 'Visible', 'off');
%     
%     gen_ecg_handle = tmp_handle;
% end

%% MECG Preparation
if choice(4)
    el = 1:min([CH_CANC, NB_EL]);
    tmp_handle = struct;
    tmp_handle.title = 'Maternal ECG cancellation';
    tmp_handle.plots = -1*ones(length(el),1);
    %set(tmp_handle.plots, 'Visible', 'off');
    %noise_handles{end+1} = tmp_handle;
        
        
    fs = 1000;
    param = out.param;
    
    for ch=1:length(el)
        peaks = adjust_mqrs_location(out.mixture(el(ch),:),out.mqrs,param.fs,0);
        peaks_1to15 = peaks(peaks<tm_idx);
        [residual, ~] = mecg_cancellation(peaks,out.mixture(el(ch),:),'TS-CERUTTI',20,2,1000,11);
        tm = 1/fs:1/fs:length(residual)/fs;        
        ecg = out.mixture(el(ch),:);
        %% MECG Plot: Maternal ECG cancellation
        %tmp_handle.plots(ch) = figure('name',sprintf('Maternal ECG cancellation, Channel %d',channels(ch)));
        tmp_handle.plots(ch) = figure('name',sprintf('Channel %d',el(ch)));
        set(tmp_handle.plots(ch), 'Visible', 'off');

        plot(tm(1:tm_idx),ecg(1:tm_idx),'LineWidth',3);
        hold on, plot(tm(1:tm_idx),ecg(1:tm_idx)-residual(1:tm_idx),'--k','LineWidth',3);
        hold on, plot(tm(1:tm_idx),residual(1:tm_idx)-1.5,'--r','LineWidth',3);
        hold on, plot(tm(peaks_1to15),ecg(peaks_1to15),'+r','LineWidth',2);
        hold off
        legend('mixture','template','residual','MQRS'); 
        title('Template subtraction for extracting the FECG');
        xlabel('Time [sec]'); ylabel('Amplitude [NU]')
        xlim([1 4]);
    end
    
    mecg_handle = tmp_handle;
    
end




end


