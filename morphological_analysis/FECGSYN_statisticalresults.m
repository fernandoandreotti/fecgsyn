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


% Alternative plots with medians for each case
count1=1;
figure
hold on
statsall = zeros(5,7,8);
symlist = {'+' 'o' '*' '.' 'x' 's' 'd'};
snrlist = [0.690196078431373 0.831372549019608 0.674509803921569;...
    0.462745098039216 0.709803921568628 0.541176470588235;...
    0.243137254901961 0.615686274509804 0.447058823529412;...
    0.0352941176470588 0.541176470588235 0.360784313725490;...
    0 0.478431372549020 0.278431372549020];


for met = {'JADEICA' 'PCA' 'tsc' 'tspca' 'tsekf' 'alms' 'arls' 'aesn'}
    eval(['stat = stats.' met{:} ';']);
    count2 = 1;
    for snr = {'snr00' 'snr03' 'snr06' 'snr09' 'snr12'}
        snrloop = eval(snr{:});
        statsall(count2,1:7,count1) = median(100*[stat(base&snrloop,1) stat(c0&snrloop,1)...
            stat(c1&snrloop,1) stat(c2&snrloop,1) stat(c3&snrloop,1) ...
            stat(c4&snrloop,1) stat(c5&snrloop,1)]);
        count2 = count2 + 1;
    end
    count1 = count1+1;
end

% plotting
count3 = linspace(-0.3,0.3,7);
for i1 = 1:size(statsall,1)
    for i2 = 1:size(statsall,2)
        plot([1:8]+ones(1,8).*count3(i2),reshape(statsall(i1,i2,:),1,8),symlist{i2},'Color',snrlist(i1,:),'MarkerSize',2*i1+3);
    end
end
LWIDTH = 1.5;
FSIZE = 15;
h=gca;
set(h, 'LineWidth',LWIDTH)
ylabel('F_1 (%)','FontSize',FSIZE)
set(gca,'XTick',[1:8])  % This automatically sets
set(gca,'XTickLabel',{'BSSica';'BSSpca';'TSc';'TSpca';'TSekf';'Alms';'Arls';'Aesn'})
set(gca,'FontSize',FSIZE);
set(findall(gcf,'-property','FontSize'),'FontSize',FSIZE);
hold off

%% Statistical tests
%F1
count = 1;
for str1 = {'base','c0','c1','c2','c3','c4','c5'}
    statsf1(count,:) = eval(['median(stats_f1(' str1{:} ',:))']);
    count = count +1;
end
p = friedman(statsf1,1,'on');
if p < 0.05
    sprintf('Group central tendency (median) are significantly different at p = %d',p)
else
    sprintf('Group central tendency (median) is NOT significantly different (p = %d)',p)
end
psig = eye(8);
hsig = eye(8);
for i = 1:8
    for j = 1:8
        [psig(i,j),hsig(i,j)] = signtest(statsf1(:,i),statsf1(:,j));
    end
end
disp(psig)
%MAE
count = 1;
for str1 = {'base','c0','c1','c2','c3','c4','c5'}
    statsmae(count,:) = eval(['median(stats_MAE(' str1{:} ',:))']);
    count = count +1;
end
p = friedman(statsmae,1,'on');
if p < 0.05
    sprintf('Group central tendency (median) are significantly different at p = %d',p)
else
    sprintf('Group central tendency (median) is NOT significantly different (p = %d)',p)
end
psig = eye(8);
hsig = eye(8);
for i = 1:8
    for j = 1:8
        [psig(i,j),hsig(i,j)] = signtest(statsmae(:,i),statsmae(:,j));
    end
end
disp(psig)