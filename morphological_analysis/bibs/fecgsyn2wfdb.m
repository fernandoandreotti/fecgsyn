% Convert all .mat files in directory to WFDB format
%
% This script converts recording from FECGSYN Matlab/Octave format to
% Physionet's WFDB format.
%
%
% NI-FECG simulator toolbox, version 1.0, February 2014
% Released under the GNU General Public License
%
% Copyright (C) 2014  Joachim Behar & Fernando Andreotti
% Oxford university, Intelligent Patient Monitoring Group - Oxford 2014
% joachim.behar@eng.ox.ac.uk, fernando.andreotti@mailbox.tu-dresden.de
%
% Last updated : 10-03-2016
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free So                                                           it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.


function fecgsyn2wfdb(lpath)
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
    signame = regexprep(cellstr([repmat('ch',size(out.mecg,1),1),num2str([1:size(out.mecg,1)]')]),'[^\w'']','');
    % save out.param as comment
    out.param = rmfield(out.param,'elpos'); % not necessary
    sn=fieldnames(out.param); % check struct
    sc=struct2cell(out.param);
    convt=cellfun(@(x) iscell(x), sc); % some cells to strings
    sa = sc(convt); idx = find(convt);
    for i = 1:length(sa)
        value = sa{i};
        switch length(value)
            case 0
                sc(idx(i)) = {'none'};
            case 1
                if isnumeric(value{:})
                    sc(idx(i)) = {value{:}};
                else
                    sc(idx(i)) = cellstr(strjoin(value));
                end
            otherwise
                if isnumeric(value{:,1})
                    sc(idx(i)) = cellstr(strjoin(cellfun(@(x) sprintf('%d ',x),value,'UniformOutput',0)));
                else
                    sc(idx(i)) = cellstr(strjoin(value));
                end
        end
        
    end
    convt=cellfun(@(x) isnumeric(x), sc); % some number need to be converted to strings
    sc(convt)=cellfun(@(x) sprintf('%d',x),sc(convt),'UniformOutput',0);
    info=cellstr(strcat(char(sn),repmat(':',size(sc)),char(sc))); % set information
    
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
        xbit=mat2wfdb(sigres,filename,fsnew,bit,[],info,gain,signame,[],false); % save in wfdb format
        ext = regexprep(signals{:},'[^a-zA-Z0-9]',''); % extension
        movefile([filename '.dat'], [filename '.' ext]); % renaming extension
    end
    flist = cellstr(arrayfun(@(x) strcat(sprintf('fqrs{%d}\n',x)),1:length(out.fqrs),'UniformOutput',0)); % list of fecg sources
    for ann = {'mqrs' flist{:}}
        annot = eval(['out.' ann{:}])/(fs/fsnew);
        ext = regexprep(ann{:},'[^a-zA-Z0-9]',''); % extension
        wrann(filename,ext,annot',repmat('N',length(annot),1));
    end
end


