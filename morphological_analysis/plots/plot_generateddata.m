cd('/media/andreotti/FetalEKG/2014.10_fecgsyn_simulations(5.0)/')
fls = dir('*.mat');     % looking for .mat (creating index)
fls =  arrayfun(@(x)x.name,fls,'UniformOutput',false);
fs_new = 250;

for i = 617:623%length(fls)
        tic
        disp(['Extracting file ' fls{i} '..'])
        % = loading data
        load(fls{i})
        disp(num2str(i))
        if isempty(out.noise)
            noise = zeros(size(out.mecg));
        else
            noise = sum(cat(3,out.noise{:}),3);
        end
        fs = out.param.fs;
        INTERV = round(0.05*fs_new);    % BxB acceptance interval
        TH = 0.3;                   % detector threshold
        REFRAC = .15;               % detector refractory period (in s)
        mixture = double(out.mecg) + sum(cat(3,out.fecg{:}),3) ...
            + noise;     % re-creating abdominal mixture
        mixture = mixture./3000;    % removing gain given during int conversion
        out.fqrs{1} = round(out.fqrs{1}/(fs/fs_new));
        out.mqrs = round(out.mqrs/(fs/fs_new));
        mixture=resample(mixture(11,:),fs_new,fs);
        plot(mixture,'b')
        xlabel('Time(s)','FontSize',14,'FontWeight','bold'), ylabel('Amplitude (mV)','FontSize',14,'FontWeight','bold')
        set(gca,'XTick',1:fs_new*10:size(mixture,2)) % This automatically sets
        set(gca,'XTickLabel',num2cell(0:10:size(mixture,2)/fs_new))
        set(gca,'FontSize',12)
        hold on
        plot(out.mqrs,mixture(out.mqrs),'dg')
        plot(out.fqrs{1},mixture(out.fqrs{1}),'or')
        legend('mixture','MQRS','FQRS')
        hold off
%         counter = 1;

        %         for k = 1:1250:length(mixture)
%             xlim([k k+1250])
% %             print('-depsc',['example_' num2str(i) '_' num2str(counter)])
%             print('-dpng',['example_' num2str(i) '_' num2str(counter)])
%             counter = counter +1 ;
%         end
%         close all
end