function FECGSYN_exp3results
% this script generates the plots for experiment 3 of the paper
%
% NI-FECG simulator toolbox, version 1.0, February 2014
% Released under the GNU General Public License
%
% Copyright (C) 2014 Joachim Behar & Fernando Andreotti
% Oxford university, Intelligent Patient Monitoring Group - Oxford 2014
% joachim.behar@eng.ox.ac.uk, fernando.andreotti@mailbox.tu-dresden.de
%
% Last updated : 30-05-2014
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License


%% Generating scatter FQT / FTQRS plots
load('/mnt/Data/Andreotti/PhD/Publications/Periodicals/2015.03 Physiol Meas - ICA breaks down/Simulator Paper/Final_Results/tqrsresult.mat')
bas = cellfun(@(x) isempty(regexp(x,'_c[0-7]','match')),fls_orig); % find out which case it depicts  baselines
% FQT
colm = [1,2;3,4];
for row = 1:2;
    count = 1;
    figure
for met = {'JADEICA' 'PCA' 'tsc' 'tspca' 'tsekf' 'alms' 'arls' 'aesn' }
    disp(met{:})
    morphtmp = morphall.(met{:})(bas,:);
    fqt = morphtmp(:,colm(row,:));
    nans.(met{:}) = sum(cellfun(@(x) sum(sum(isnan(cell2mat(x)))),fqt));
    tot.(met{:}) = sum(cellfun(@(x) numel(cell2mat(x)),fqt));
    fqtmed = cellfun(@(x) abs(cell2mat(x)),fqt,'UniformOutput',0);
    fqtmed = cellfun(@(x) nanmedian(x),fqtmed,'UniformOutput',0);
    fqtmed = [cell2mat(fqtmed(:,1)')' cell2mat(fqtmed(:,2)')'];
    fqtmed(any(isnan(fqtmed), 2),:)=[];
    [h,coe]=ttest(fqtmed(:,2),fqtmed(:,1));
    p.(met{:}) = coe;
    srcnan = nans.(met{:})(1); srctot = tot.(met{:})(1);
    timenan = nans.(met{:})(2); timetot = tot.(met{:})(2);
    if h
        fprintf('Null hypothesis can be rejected with p= %d \n',p.(met{:}));
    else
        fprintf('Null hypothesis CANNOT be rejected');
    end
    fprintf('NaNs on source domain:  %3.2f percent \n',srcnan/srctot*100);
    fprintf('NaNs on time domain:  %3.2f percent \n',timenan/timetot*100);
    % plot
    subplot(2,4,count)
    x = fqtmed(:,2); y = fqtmed(:,1);
    scatter(x,y,12,'filled')
     % regression
    a = polyfit(x,y,1);
    yfit = polyval(a,x);
    yresid = y - yfit;
    SSresid = sum(yresid.^2);
    SStotal = (length(y)-1) * var(y);
    rsq = 1 - SSresid/SStotal;
    hold on; plot(x,yfit,'r-'); 
    if row ==1
    text(220,140,sprintf('r^2 = %2.4f',rsq))
    xlim([120,260]),ylim([120,260])
    xlabel('FQT reference (ms)')
    ylabel('FQT test (ms)')
    else
       text(1,0.2,sprintf('r^2 = %2.4f',rsq))
        xlim([0 1.5]),ylim([0 5])
        xlabel('T/R reference (ms)')
        ylabel('T/R test (ms)')
    end
    
    title(met{:})

    count = count +1;
end
end
%= Distribution of ICA
% FQT intervals
exp3dist = exp3dist1;

    
    exp3med = cellfun(@(x) nanmax(cell2mat(x)),exp3dist,'UniformOutput',0);
    exp3med = [cell2mat(exp3med(:,1)')' cell2mat(exp3med(:,2)')'];
    [h,p]=ttest(exp3med(:,2),exp3med(:,1));
   
    x = exp3med(:,2); y = exp3med(:,1);
    scatter(x,y,12,'filled')
   
    
    % FTh
    exp3dist = exp3dist1;
    nans = sum(cellfun(@(x) sum(sum(isnan(cell2mat(x)))),exp3dist));
    tot = sum(cellfun(@(x) numel(cell2mat(x)),exp3dist));
    srcnan = nans(1); srctot = tot(1);
    timenan = nans(2); timetot = tot(2);
    exp3med = cellfun(@(x) nanmedian(abs(cell2mat(x))),exp3dist,'UniformOutput',0);
    exp3med = [cell2mat(exp3med(:,3)')' cell2mat(exp3med(:,4)')'];
    [h,p]=ttest(exp3med(:,2),exp3med(:,1));
    if h
        fprintf('Null hypothesis can be rejected with p= %d \n',p);
    else
        fprintf('Null hypothesis CANNOT be rejected');
    end
    fprintf('NaNs on source domain:  %3.2f percent \n',srcnan/srctot*100);
    fprintf('NaNs on time domain:  %3.2f percent \n',timenan/timetot*100);
     x = exp3med(:,2); y = exp3med(:,1);
    scatter(x,y,12,'filled')
    % regression
    A = [x,y];
    A(any(isnan(A), 2),:)=[];
    x = A(:,1); y = A(:,2);
    p = polyfit(x,y,1);
    yfit = polyval(p,x);
    yresid = y - yfit;
    SSresid = sum(yresid.^2);
    SStotal = (length(y)-1) * var(y);
    rsq = 1 - SSresid/SStotal;
    hold on; plot(x,yfit,'r-'); 
    text(1.5,0.1,sprintf('r^2 = %2.4f',rsq))
   % xlim([120,260]),ylim([120,260])
    xlabel('FT_h in time domain (n.u.)')
    ylabel('FT_h in source domain (n.u.)')
    
    % Case by case methods against each other
    % Generate Table
    res = struct('qt',[],'th',[]);
    qt = []; th = [];
    % FQT
    for met = {'JADEICA' 'PCA' 'tsc' 'tspca' 'tsekf' 'alms' 'arls' 'aesn' }
        tmp = morphall.(met{:});
        res.qt=cell(1750,1);
        for i = 1:1750
            res.qt{i} = cell2mat(tmp{i,1})-cell2mat(tmp{i,2});
        end
        %         stat = bsxfun(@minus,morphall.(met{:})(:,col),morphall.(met{:})(:,col+1));
        res.qt = cellfun(@(x) median(nanmin(x)),res.qt);
        res.qtstd = cellfun(@(x) std(nanmin(x)),res.qt);
        qt = [qt nanmedian(res.qt)];
        qtstd = [qtstd nanmedian(res.qtstd)];
    end
    
    % FTh
    for met = {'JADEICA' 'PCA' 'tsc' 'tspca' 'tsekf' 'alms' 'arls' 'aesn' }
        tmp = morphall.(met{:});
        res.th=cell(1750,1);
        for i = 1:1750
            res.th{i} = cell2mat(tmp{i,3})./cell2mat(tmp{i,4});
        end
        res.th = cellfun(@(x) nanmedian(nanmin(x-1)),res.th);
        th = [th nanmedian(res.th)];
    end
    
    