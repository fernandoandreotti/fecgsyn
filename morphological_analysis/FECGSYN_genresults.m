function [outsig,qrsmethod] = FECGSYN_genresults(path_orig,path_ext,fs,ch,debug)
%
% Input:
%  path_orig:       Path for original dataset
%  path_ext:        Path for extracted dataset
%  fs:              Sampling frequency
% this script generates a series of abdominal mixtures, containing i) a
% stationary case and ii) non-stationary case (when adding breathing
% effects, foetal movement etc.).
%
%
% NI-FECG simulator toolbox, version 1.0, February 2014
% Released under the GNU General Public License
%
% Copyright (C) 2014  Joachim Behar & Fernando Andreotti
% Oxford university, Intelligent Patient Monitoring Group - Oxford 2014
% joachim.behar@eng.ox.ac.uk, fernando.andreotti@mailbox.tu-dresden.de
%
% Last updated : 09-08-2014
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

%% == Parameters
INTERV = round(0.05*fs);    % BxB acceptance interval
TEMP_SEC = round(60*fs);    % samples used for building templates
morph = 1;                  % turn on/off morphological analysis
%% Run through extracted datasets
cd(path_orig)
slashchar = char('/'*isunix + '\'*(~isunix));
fls_orig = dir('*.mat');     % looking for .mat (creating index)
fls_orig =  arrayfun(@(x)x.name,fls_orig,'UniformOutput',false);
cd(path_ext)
fls_ext = dir('*.mat');     % looking for .mat (creating index)
fls_ext =  arrayfun(@(x)x.name,fls_ext,'UniformOutput',false);

%
stats_ica = zeros(length(fls_orig),4);
stats_pca = zeros(length(fls_orig),4);
stats_tsc = zeros(length(fls_orig),4);
stats_tspca = zeros(length(fls_orig),4);
stats_tsekf = zeros(length(fls_orig),4);
stats_alms = zeros(length(fls_orig),4);
stats_arls = zeros(length(fls_orig),4);
stats_aesn = zeros(length(fls_orig),4);

morph_ica = zeros(length(fls_orig),2);
morph_tspca = zeros(length(fls_orig),2);
morph_tskf = zeros(length(fls_orig),2);
morph_aesn = zeros(length(fls_orig),2);


for i = 1:length(fls_ext)
    disp(fls_ext{i})
    [rec,met] = strtok(fls_ext(i),'_');
    file = strcat(path_ext,fls_ext(i)); 
    load(file{:})
    origrec = str2double(rec{:}(4:end));
    file = strcat(path_orig,fls_orig(origrec)); 
    load(file{:});
    [elif,~]=strtok(file{:}(end:-1:1),slashchar);
    disp(elif(end:-1:1))
    switch met{:}(2:end-4)
        case 'JADEICA'
            % generating statistics
            [F1,MAD,PPV,SE] = Bxb_compare(out.fqrs{1},fqrs,INTERV);
            stats_ica(origrec,:) = [F1,MAD,PPV,SE];
            %             if morph
            %             % generating reference template
            %             % begin(workaround) while not regenerating data at 250 Hz
            %             for j = 1:length(ch)
            %                 fecg(j,:) = double(out.fecg{1}(ch(j),:));
            %             end
            %             % end(workaround)
            %             %  fecg = double(out.fecg{1}(ch,:));   when workaround is undone
            %             W = fecg*pinv(outdata);
            %             srcfecg = W*fecg;
            %
            %             % templates generation
            %             for j = 1:length(maxch)
            %                 ssamp = (j-1)*round(length(outdata)/length(maxch))+1;
            %                 endsamp = j*round(length(outdata)/length(maxch))-1;
            %                 qrstmp = fqrs(fqrs>ssamp&fqrs<endsamp);
            %
            %             end
            %             end
            clear fqrs F1 MAD PPV SE
        case 'PCA'
            % generating statistics
            [F1,MAD,PPV,SE] = Bxb_compare(out.fqrs{1},fqrs,INTERV);
            stats_pca(origrec,:) = [F1,MAD,PPV,SE];
            clear fqrs F1 MAD PPV SE
        case 'tsc'
            % discarding channels that are not the best
            fqrs = fqrs{maxch};
            residual = double(residual(maxch,:))./3; %in mV
            fecg = double(out.fecg{1}(ch(maxch),:))./3; %in mV
            % generating QRS detection statistics
            [F1,MAD,PPV,SE] = Bxb_compare(out.fqrs{1},fqrs,INTERV);
            stats_tsc(origrec,:) = [F1,MAD,PPV,SE];
                        
            clear fqrs F1 MAD PPV SE
            
        case 'tspca'
            % discarding channels that are not the best
            fqrs = fqrs{maxch};
            residual = residual(maxch,:);
            fecg = double(out.fecg{1}(maxch,:));
            % generating statistics
            [F1,MAD,PPV,SE] = Bxb_compare(out.fqrs{1},fqrs,INTERV);
            stats_tspca(origrec,:) = [F1,MAD,PPV,SE];
            
            % morphological statistics
            if morph
                qt_err = []; theight_err = [];
                for j = 1:TEMP_SEC:length(residual)
                    % checking borders
                    if j+TEMP_SEC > length(residual)
                        endsamp = length(residual);
                    else
                        endsamp = j + TEMP_SEC;
                    end
                    % qrs complexes in interval
                    qrstmp = fqrs(fqrs>j&fqrs<endsamp)-j;
                    % abdominal signal template
                    temp_abdm = FECGSYN_tgen(residual(j:endsamp),qrstmp,debug);
                    % reference template
                    temp_ref = FECGSYN_tgen(fecg(j:endsamp),qrstmp,debug);
                    temp_abdm = temp_abdm.avg; temp_ref = temp_ref.avg;
                    % evaluating morphological features
                    [qt_err(end+1),theight_err(end+1)] = FECGSYN_manalysis(temp_abdm,temp_ref,fs);
                end
            end
            morph_tspca(end+1,:) = [mean(qt_err); mean(theight_err)];
            
            clear fecg residual fqrs F1 MAD PPV SE qt_err theight_err
        case 'tsekf'
            % discarding channels that are not the best
            fqrs = fqrs{maxch};
            residual = residual(maxch,:);
            fecg = double(out.fecg{1}(maxch,:));
            % generating statistics
            [F1,MAD,PPV,SE] = Bxb_compare(out.fqrs{1},fqrs,INTERV);
            stats_tsekf(origrec,:) = [F1,MAD,PPV,SE];
            
            % morphological statistics
            if morph
                qt_err = []; theight_err = [];
                for j = 1:TEMP_SEC:length(residual)
                    % checking borders
                    if j+TEMP_SEC > length(residual)
                        endsamp = length(residual);
                    else
                        endsamp = j + TEMP_SEC;
                    end
                    % qrs complexes in interval
                    qrstmp = fqrs(fqrs>j&fqrs<endsamp)-j;
                    % abdominal signal template
                    temp_abdm = FECGSYN_tgen(residual(j:endsamp),qrstmp,debug);
                    % reference template
                    temp_ref = FECGSYN_tgen(fecg(j:endsamp),qrstmp,debug);
                    temp_abdm = temp_abdm.avg; temp_ref = temp_ref.avg;
                    % evaluating morphological features
                    [qt_err(end+1),theight_err(end+1)] = FECGSYN_manalysis(temp_abdm,temp_ref,fs);
                end
            end
            morph_tskf(end+1,:) = [mean(qt_err); mean(theight_err)];
            
            clear fecg residual fqrs F1 MAD PPV SE qt_err theight_err
        case 'alms'
            % discarding channels that are not the best
            fqrs = fqrs{maxch};
            residual = residual(maxch,:);
            fecg = double(out.fecg{1}(maxch,:));
            % generating statistics
            [F1,MAD,PPV,SE] = Bxb_compare(out.fqrs{1},fqrs,INTERV);
            stats_alms(origrec,:) = [F1,MAD,PPV,SE];
            clear fecg residual fqrs F1 MAD PPV SE
        case 'arls'
            % discarding channels that are not the best
            fqrs = fqrs{maxch};
            residual = residual(maxch,:);
            fecg = double(out.fecg{1}(maxch,:));
            % generating statistics
            [F1,MAD,PPV,SE] = Bxb_compare(out.fqrs{1},fqrs,INTERV);           
            stats_arls(origrec,:) = [F1,MAD,PPV,SE];
             clear fecg residual fqrs F1 MAD PPV SE
        case 'aesn'
            % discarding channels that are not the best
            fqrs = fqrs{maxch};
            residual = residual(maxch,:);
            fecg = double(out.fecg{1}(maxch,:));
            % generating statistics
            [F1,MAD,PPV,SE] = Bxb_compare(out.fqrs{1},fqrs,INTERV);
            stats_aesn(origrec,:) = [F1,MAD,PPV,SE];
            
            % morphological statistics
            if morph
                qt_err = []; theight_err = [];
                for j = 1:TEMP_SEC:length(residual)
                    % checking borders
                    if j+TEMP_SEC > length(residual)
                        endsamp = length(residual);
                    else
                        endsamp = j + TEMP_SEC;
                    end
                    % qrs complexes in interval
                    qrstmp = fqrs(fqrs>j&fqrs<endsamp)-j;
                    % abdominal signal template
                    temp_abdm = FECGSYN_tgen(residual(j:endsamp),qrstmp,debug);
                    % reference template
                    temp_ref = FECGSYN_tgen(fecg(j:endsamp),qrstmp,debug);
                    temp_abdm = temp_abdm.avg*1000; % looks like adaptive filters are scalling them
                    temp_ref = temp_ref.avg;
                    % evaluating morphological features
                    [qt_err(end+1),theight_err(end+1)] = FECGSYN_manalysis(temp_abdm,temp_ref,fs);
                end
            end
            morph_aesn(end+1,:) = [mean(qt_err); mean(theight_err)];
            
            clear fecg residual fqrs F1 MAD PPV SE qt_err theight_err
    end
end

%% Statistics Generation
LWIDTH = 1.5;
FSIZE = 15;
% F1
stats_f1 = 100*[stats_ica(:,1) stats_pca(:,1) stats_tsc(:,1) stats_tspca(:,1) ...
    stats_tsekf(:,1) stats_alms(:,1) stats_arls(:,1) stats_aesn(:,1)];
h = boxplot(stats_f1,{[repmat({'BSSica'},1,length(fls_orig)) repmat({'BSSpca'},1,length(fls_orig)) ...
    repmat({'TSc'},1,length(fls_orig)) repmat({'TSpca'},1,length(fls_orig)) ...
    repmat({'TSekf'},1,length(fls_orig)) repmat({'Alms'},1,length(fls_orig)) ...
    repmat({'Arls'},1,length(fls_orig)) repmat({'Aesn'},1,length(fls_orig))]});
set(h, 'LineWidth',LWIDTH)
ylabel('F1 (%)','FontSize',FSIZE)
% MAD
stats_MAD = [stats_ica(:,2) stats_pca(:,2) stats_tsc(:,2) stats_tspca(:,2) ...
    stats_tsekf(:,2) stats_alms(:,2) stats_arls(:,2) stats_aesn(:,2)];
h = boxplot(stats_MAD,{[repmat({'ICA'},1,length(fls_orig)) repmat({'PCA'},1,length(fls_orig)) ...
    repmat({'TSc'},1,length(fls_orig)) repmat({'TSpca'},1,length(fls_orig)) ...
    repmat({'EKF'},1,length(fls_orig)) repmat({'LMS'},1,length(fls_orig)) ...
    repmat({'RLS'},1,length(fls_orig)) repmat({'ESN'},1,length(fls_orig))]});
set(h, 'LineWidth',LWIDTH)
ylabel('MAD (ms)','FontSize',FSIZE)


% Plot about cases

% This script plots boxplots with 2 groups
% Boxplot multicolor
% % stats_ica(:,[1 3 4]) = 100.*stats_ica(:,[1 3 4]);
% % stats_tsc(:,[1 3 4]) = 100.*stats_tsc(:,[1 3 4]);
% % LINE_WIDTH = 1.2;
% % FONT_SIZE = 12;
% % MARKER_SIZE = 7;
% % fig1=figure(1)
% % base = cellfun(@(x) ~isempty(regexp(x, '_l\d.mat$', 'match')), fls_orig);
% % c0 = cellfun(@(x) ~isempty(regexp(x, '_c0.mat$', 'match')), fls_orig);
% % c1 = cellfun(@(x) ~isempty(regexp(x, '_c1.mat$', 'match')), fls_orig);
% % c2 = cellfun(@(x) ~isempty(regexp(x, '_c2.mat$', 'match')), fls_orig);
% % c3 = cellfun(@(x) ~isempty(regexp(x, '_c3.mat$', 'match')), fls_orig);
% % c4 = cellfun(@(x) ~isempty(regexp(x, '_c4.mat$', 'match')), fls_orig);
% % c5 = cellfun(@(x) ~isempty(regexp(x, '_c5.mat$', 'match')), fls_orig);
% % 
% % stats_ica = stats_ica(:,1);
% % stats_tsc = stats_tsc(:,1);
% % 
% % N = sum(base)+sum(c0)+sum(c1)+sum(c2)+sum(c3)+sum(c4)+sum(c5);
% % bh=boxplot([stats_ica(base); stats_ica(c0); stats_ica(c1); stats_ica(c2); stats_ica(c3); stats_ica(c4); stats_ica(c5);...
% %     stats_tsc(base); stats_ica(c0); stats_tsc(c1); stats_tsc(c2); stats_tsc(c3); stats_tsc(c4); stats_tsc(c5)], ...
% %     {[repmat({'Base'},1,sum(base)) repmat({'Case 0'},1,sum(c0)) repmat({'Case1'},1,sum(c1)) repmat({'Case2'},1,sum(c2)) repmat({'Case3'},1,sum(c3)) repmat({'Case4'},1,sum(c4)) repmat({'Case5'},1,sum(c5)) ...
% %     repmat({'Base'},1,sum(base)) repmat({'Case 0'},1,sum(c0)) repmat({'Case1'},1,sum(c1)) repmat({'Case2'},1,sum(c2)) repmat({'Case3'},1,sum(c3)) repmat({'Case4'},1,sum(c4)) repmat({'Case5'},1,sum(c5)) ] ...
% %     [repmat({'ICA'},N,1);repmat({'TS'},N,1)]},'colors',repmat('rk',1,6),'factorgap',10,'labelverbosity','minor','labelorientation','inline')
% % set(bh,'linewidth',LINE_WIDTH);
% % ylabel('F1 (in %)','FontSize',FONT_SIZE)
% % xlabel('Recording','FontSize',FONT_SIZE)
% % set(gca,'FontSize',FONT_SIZE)
% % set(findall(gcf,'type','text'),'fontSize',FONT_SIZE);
% % %     ylim([0 105])
% % xl=xlabel('Recording Number','FontSize',FONT_SIZE);
% % pos=get(xl,'Pos');
% % set(xl,'Pos',[pos(1) pos(2)-30 pos(3)])
% % save2pdf('boxplot',fig1,600)
% % 
% % % Plot about SNRmn
% % a = cellfun(@(x) strsplit(x,'_snr'), fls_orig,'UniformOutput',0);
% % 
% % b = cellfun(@(x) length(x)>1,a);
% % 


