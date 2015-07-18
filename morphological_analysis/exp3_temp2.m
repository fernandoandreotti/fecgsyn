% This script runs through comparison between QT and Th values for applying
% or not the 
cd('/home/andreotti/Desktop/Simulator Paper/generalresults')
fls = dir('*.mat');     % looking for .mat (creating index)
fls =  arrayfun(@(x)x.name,fls,'UniformOutput',false);
morphall = struct('JADEICA',[],'PCA',[],'tsc',[],'tspca',[],'tsekf',[],'alms',[],'arls',[],'aesn',[]);
for i = 1:length(fls)
    load(fls{i})
    
end