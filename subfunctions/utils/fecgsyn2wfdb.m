function fecgsyn2wfdb(lpath)
% function fecgsyn2wfdb(argument)
% Convert all .mat files in directory to WFDB format
%
% This script converts recording(s) from fecgsyn's Matlab/Octave format to
% Physionet's WFDB format.
%
%  Input:
%  argument         either a local path (string) or the "out" structure
%                   intern from fecgsyn
%
%
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
%
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

slashchar = char('/'*isunix + '\'*(~isunix));

outpath = [lpath slashchar 'wfdb' slashchar];
fls = dir([lpath  '*.mat']);   % looking for .mat (creating index)
fls =  arrayfun(@(x)x.name,fls,'UniformOutput',false);
gain = 3000; % files are converted to integers with a 3000 gain
bit = 16;
fsnew = 250; % new sampling rate
cd(outpath)

for i = 1:length(fls)
    load([lpath fls{i}],'out') % load dataset
    % getting information
    filename = ['sub' fls{i}(8:end-4)];
    fs = out.param.fs;
    signame = cellstr(regexprep(cellstr([repmat('ch',size(out.mecg,1),1),...
        num2str([1:size(out.mecg,1)]')]),'[^\w'']',''));
    % save out.param as comment
    out.param = rmfield(out.param,'elpos'); % not necessary
    out.param = rmfield(out.param,'fs');    % not necessary, specified in wfdb header
    if isfield(out.param,'noise_fct'), out.param = ...
            rmfield(out.param,'noise_fct');end % applied noise function can be derived from signal itself
    sn=fieldnames(out.param); % check struct
    sc=struct2cell(out.param);
    convt=cellfun(@(x) iscell(x), sc); % some cells to strings
    % The struct is fairly heterogeneous and requires this kind of
    % workarounds to convert every data type possible and human readable
    % strings
    sa = sc(convt); idx = find(convt);
    for i = 1:length(sa)
        value = sa{i};
        switch length(value)
            case 0
                sc(idx(i)) = {'none'};
            case 1
                if isnumeric(value{:})
                    sc(idx(i)) = {sprintf('%2.4f ',value{:})};
                else
                    sc(idx(i)) = cellstr(strjoin(value));
                end
            otherwise
                if isnumeric(value{:,1})
                    sc(idx(i)) = cellstr(strjoin(cellfun(@(x) sprintf('%2.4f ',x),value,'UniformOutput',0)));
                else
                    sc(idx(i)) = cellstr(strjoin(value));
                end
        end
        
    end
    convt=cellfun(@(x) isnumeric(x), sc); % some number need to be converted to strings
    % converting remaining floats/integers
    convt = find(convt);
    for idx = 1:length(convt)        
        if rem(sc{convt(idx)},1)==0
            sc(convt(idx))= cellstr(sprintf('%d ',sc{convt(idx)}));
        else
            sc(convt(idx))= cellstr(sprintf('%2.4f ',sc{convt(idx)}));
        end
    end
    info=cellstr(strcat(char(sn),repmat(':',size(sc)),char(sc))); % set information
    info(end+1,1) = cellstr(sprintf('nfetus:%d',length(out.fecg)));
    info(end+1,1) = cellstr(sprintf('nnoise:%d',length(out.noise)));
    % converting data
    nlist = (arrayfun(@(x) strcat(sprintf('noise{%d}\n',x)),1:length(out.noise),'UniformOutput',0)); % list of noise sources
    flist = cellstr(arrayfun(@(x) strcat(sprintf('fecg{%d}\n',x)),1:length(out.fecg),'UniformOutput',0)); % list of fecg sources
    
    for signals = {'mecg' flist{:} nlist{:}}
        % resampling
        sig = eval(['out.' signals{:}])';
        sigres = zeros(size(sig,1)/(fs/fsnew),size(sig,2));
        for ch = 1:size(sigres,2)
            sigres(:,ch) = resample(double(sig(:,ch)),fsnew,fs);
        end
        ext = regexprep(signals{:},'[^a-zA-Z0-9]',''); % extension
%        wrsamp(tm2,ref_sig',recordName,FS_ECGPU,gain,'')
    xbit=mat2wfdb(sigres,[filename '_' ext],fsnew,bit,{'nu'},info,gain,signame,0,false); % save in wfdb format
    end
    flist = cellstr(arrayfun(@(x) strcat(sprintf('fqrs{%d}\n',x)),1:length(out.fqrs),'UniformOutput',0)); % list of fecg sources
    for ann = {'mqrs' flist{:}}
        annot = round(eval(['out.' ann{:}])/(fs/fsnew));
        ext = regexprep(ann{:},'[^a-zA-Z0-9]',''); % extension
        wrann([filename '_mecg'],ext,annot',repmat('N',length(annot),1));
    end
end


