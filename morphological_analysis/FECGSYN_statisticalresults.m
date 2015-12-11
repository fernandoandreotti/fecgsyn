% Loading data
load('D:\Users\Andreotti\Dropbox\sharelatex\[Andreotti et al 2016] Standardising FECGSYN (in review)\fecgsyn_morphology\mats\exp2.mat')
fls_orig(1751:end) = []

% Marking cases
c0 = cellfun(@(x) ~isempty(regexp(x,'.c0','ONCE')),fls_orig);
c1 = cellfun(@(x) ~isempty(regexp(x,'.c1','ONCE')),fls_orig);
c2 = cellfun(@(x) ~isempty(regexp(x,'.c2','ONCE')),fls_orig);
c3 = cellfun(@(x) ~isempty(regexp(x,'.c3','ONCE')),fls_orig);
c4 = cellfun(@(x) ~isempty(regexp(x,'.c4','ONCE')),fls_orig);
c5 = cellfun(@(x) ~isempty(regexp(x,'.c5','ONCE')),fls_orig);
base = ~(c0|c1|c2|c3|c4|c5);
snr00 = cellfun(@(x) ~isempty(regexp(x,'.snr00dB','ONCE')),fls_orig);
snr03 = cellfun(@(x) ~isempty(regexp(x,'.snr03dB','ONCE')),fls_orig);
snr06 = cellfun(@(x) ~isempty(regexp(x,'.snr06dB','ONCE')),fls_orig);
snr09 = cellfun(@(x) ~isempty(regexp(x,'.snr09dB','ONCE')),fls_orig);
snr12 = cellfun(@(x) ~isempty(regexp(x,'.snr12dB','ONCE')),fls_orig);

symlist = {'o' '^' '<' 'v' 's' '>' 'd' }; % symbol selection for plot
snrlist = zeros(7,5,3);
% color selection for plot, different colors for each case
snrlist(1,:,:) = [[0.800000000000000 0.800000000000000 0.800000000000000] ;...
   [0.600000000000000 0.600000000000000 0.600000000000000];...
   [0.400000000000000 0.400000000000000 0.400000000000000];...
   [0.200000000000000 0.200000000000000 0.200000000000000];...
   0 0 0]; % color selection for plot
snrlist(2,:,:) = [0.784313725490196 0.635294117647059 0.784313725490196;...
   0.690196078431373 0.486274509803922 0.682352941176471;...
    0.611764705882353 0.352941176470588 0.596078431372549;...
    0.572549019607843 0.282352941176471 0.556862745098039;...
    0.505882352941176 0.101960784313725 0.470588235294118]; % color selection for plot

snrlist(3,:,:) = [[0.658823529411765 0.686274509803922 0.780392156862745] ;...
   [0.466666666666667 0.498039215686275 0.619607843137255];...
   [0.317647058823529 0.364705882352941 0.501960784313726];...
   [0.184313725490196 0.250980392156863 0.403921568627451];...
   [0.043137254901961 0.164705882352941 0.317647058823529]]; % color selection for plot

snrlist(4,:,:) = [[0.984313725490196 0.913725490196078 0.745098039215686] ;...
   [0.968627450980392 0.823529411764706 0.576470588235294];...
   [0.952941176470588 0.705882352941176 0.443137254901961];...
   [0.937254901960784 0.611764705882353 0.317647058823529];...
   [0.909803921568627 0.482352941176471 0.078431372549020]]; % color selection for plot

snrlist(5,:,:) = [0.560784313725490 0.764705882352941 0.596078431372549;...
    0.462745098039216 0.709803921568628 0.541176470588235;...
    0.243137254901961 0.615686274509804 0.447058823529412;...
    0.0352941176470588 0.541176470588235 0.360784313725490;...
    0 0.478431372549020 0.278431372549020]; % color selection for plot

snrlist(6,:,:) = [[0.772549019607843 0.725490196078431 0.854901960784314] ;...
   [0.623529411764706 0.552941176470588 0.745098039215686];...
   [0.509803921568627 0.423529411764706 0.654901960784314];...
   [0.411764705882353 0.298039215686275 0.576470588235294];...
   [0.317647058823529 0.160784313725490 0.498039215686275]]; % color selection for plot

snrlist(7,:,:) = [[0.721568627450980 0.776470588235294 0.901960784313726] ;...
   [0.529411764705882 0.631372549019608 0.823529411764706];...
   [0.380392156862745 0.521568627450980 0.752941176470588];...
   [0.203921568627451 0.435294117647059 0.698039215686275];...
   [0 0.349019607843137 0.639215686274510]]; % color selection for plot

%% Alternative plots (after review) with medians for each case, instead of boxplots

% F1
clear statsf1 statsmae
statsf1 = zeros(5,7,8); statsmae = statsf1;
count1 = 1;
for met = {'JADEICA' 'PCA' 'tsc' 'tspca' 'tsekf' 'alms' 'arls' 'aesn'}
    eval(['stat = stats.' met{:} ';']);
    count2 = 1;
    for snr = {'snr00' 'snr03' 'snr06' 'snr09' 'snr12'}
        snrloop = eval(snr{:});
        statsf1(count2,1:7,count1) = median(100*[stat(base&snrloop,1) stat(c0&snrloop,1)...
            stat(c1&snrloop,1) stat(c2&snrloop,1) stat(c3&snrloop,1) ...
            stat(c4&snrloop,1) stat(c5&snrloop,1)]); % F1
         statsmae(count2,1:7,count1) = median([stat(base&snrloop,2) stat(c0&snrloop,2)...
            stat(c1&snrloop,2) stat(c2&snrloop,2) stat(c3&snrloop,2) ...
            stat(c4&snrloop,2) stat(c5&snrloop,2)]);% MAE
        count2 = count2 + 1; 
    end
    count1 = count1+1;
end

%'statsf1' has three dimensions: SNR, Cases, Methods
% plotting
figure
hold on
count3 = linspace(-0.7,0.7,7);
for i1 = size(statsf1,1):-1:1
    for i2 = 1:size(statsf1,2)
        plot([1:2:16]+ones(1,8).*count3(i2),reshape(statsf1(i1,i2,:),1,8),symlist{i2},'Color',snrlist(i2,i1,:),'MarkerSize',round(1.5*i1+3),'MarkerFaceColor',snrlist(i2,i1,:));
    end
end
legend('baseline','case 0','case 1','case 2','case 3','case 4','case 5')
LWIDTH = 1.5;
FSIZE = 15;
h=gca;
set(h, 'LineWidth',LWIDTH)
ylabel('F_1 (%)','FontSize',FSIZE)
set(gca,'XTick',[1:2:16])  % This automatically sets
set(gca,'XTickLabel',{'BSS_{ica}';'BSS_{pca}';'TS_c';'TS_{pca}';'TS_{ekf}';'AM_{lms}';'AM_{rls}';'AM_{esn}'})
set(gca,'FontSize',FSIZE);
set(findall(gcf,'-property','FontSize'),'FontSize',FSIZE);
box on
hold off

% MAE
figure
hold on
% plotting
for i1 = size(statsmae,1):-1:1
    for i2 = 1:size(statsmae,2)
        plot([1:2:16]+ones(1,8).*count3(i2),reshape(statsmae(i1,i2,:),1,8),symlist{i2},'Color',snrlist(i2,i1,:),'MarkerSize',round(1.5*i1+3),'MarkerFaceColor',snrlist(i2,i1,:));
    end
end
legend('baseline','case 0','case 1','case 2','case 3','case 4','case 5')
h=gca;
set(h, 'LineWidth',LWIDTH)
ylabel('MAE (ms)','FontSize',FSIZE)
set(gca,'XTick',[1:2:16])  % This automatically sets
set(gca,'XTickLabel',{'BSS_{ica}';'BSS_{pca}';'TS_c';'TS_{pca}';'TS_{ekf}';'AM_{lms}';'AM_{rls}';'AM_{esn}'})
set(gca,'FontSize',FSIZE);
set(findall(gcf,'-property','FontSize'),'FontSize',FSIZE);
box on
hold off
clear count1 count2 count3 i1 i2 h 
%% Statistical tests
%F1
% re-calculating to include baseline
statsf1(:,1,:) = []; statsmae(:,1,:) = []; % removing baseline since it does not make sense on significance analysis
count1 = 1;
for var = {'statsf1' 'statsmae'}
    stastuse = eval(var{:});
    figure
    for snr = 1:size(stastuse,1)
        tempstat = reshape(stastuse(snr,:,:),[],8);
        p(snr,1) = friedman(tempstat,1,'off');
        p(snr,2) = friedman(tempstat',1,'off');
        for i = 1:8
            for j = 1:8
                [psig(snr,i,j),hsig(snr,i,j)] = signtest(tempstat(:,i),tempstat(:,j));
            end
        end
        gridp = reshape(psig(snr,:,:),8,8);
        gridp(gridp >= 0.05) = 0;
        gridp(gridp < 0.01&gridp>0) =  2;
        gridp(gridp < 0.05&gridp>0) =  1;
        subplot(1,5,snr)
        pcolor([gridp NaN(8,1);NaN(1,9)])%,[0 2])
        colormap(flipud(gray))
        %xlim([0 8]),ylim([0 8])
        set(gca,'xtick', linspace(0.5,8.5,8), 'ytick', linspace(0.5,8.5,8));
        set(gca,'xticklabel', {[1:8]}, 'yticklabel', {[1:8]});
        axis square
        %     set(gca,'xgrid', 'on', 'ygrid', 'on', 'gridlinestyle', '-', 'xcolor', 'k', 'ycolor', 'k');
        grid on
    end
end

disp(psig)
%MAE
