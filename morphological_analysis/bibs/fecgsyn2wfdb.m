% Convert all .mat files in directory to WFDB format
%
% Requires WFDB-APP 

function fecgsyn2wfdb(lpath)
slashchar = char('/'*isunix + '\'*(~isunix));

outpath = [lpath slashchar 'wfdb' slashchar]
fls = dir([lpath  '*.mat']);   % looking for .mat (creating index)
fls =  arrayfun(@(x)x.name,fls,'UniformOutput',false);

for i = 1:length(fls)
wrsamp(tm2,ref_sig',recordName,FS_ECGPU,gain,'')
wrann(recordName,'qrs',qrsref',repmat('N',20,1));