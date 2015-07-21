% This script runs through comparison between QT and Th values for applying
% or not the 
% cd('/home/andreotti/Desktop/Simulator Paper/generalresults')
% fls = dir('*.mat');     % looking for .mat (creating index)
% fls =  arrayfun(@(x)x.name,fls,'UniformOutput',false);
% morphall = struct('JADEICA',[],'PCA',[],'tsc',[],'tspca',[],'tsekf',[],'alms',[],'arls',[],'aesn',[]);
% for met = {'JADEICA','PCA','tsc','tspca','tsekf','alms','arls','aesn'}
%        morphall.(met{:}) = cell(1750,7);
%  end
% for i = 1:length(fls)
%     load(fls{i})
%     for met = {'JADEICA','PCA','tsc','tspca','tsekf','alms','arls','aesn'}
%         idx = cellfun(@(x) ~isempty(x),morph.(met{:}));
%         morphall.(met{:})(idx) = morph.(met{:})(idx);
%     end
%         
% end

% Plotting
counter = 1;
statFQT = [];
statFTH = [];
group = [];

for met = {'JADEICA','PCA','tsc','tspca','tsekf','alms','arls','aesn'}
    % FQT
    statFQT(end+1:end+1750) = cell2mat(morphall.(met{:})(:,5));
    
    % FTH
    statFTH(end+1:end+1750) = cell2mat(morphall.(met{:})(:,6));
    
    group = [group repmat(counter,1,1750)];
    counter1 = counter1 + 1;
end
boxplot(statFQT,group)