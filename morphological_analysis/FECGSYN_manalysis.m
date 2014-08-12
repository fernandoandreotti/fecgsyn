function [qt,theight] = FECGSYN_manalysis(test_temp,ref_temp,fs)
% This function calculates morphological features form signals given two
% templates (reference and test). Statistics are give as %.
% 
% Input:
%  test_temp:       Template to be tested
%  path_ext:        Path for extracted dataset
%  fs:              Sampling frequency
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

% == generate signal from template (necessary for ecgpuwave)
fsnew = 2*250*round(250/fs);
test_sig = repmat(test_temp,1,20);
test_sig = test_sig/max(abs(test_sig));
test_sig = resample(test_sig,fsnew,fs);
tm = 1:length(test_sig)-1;
wrsamp(tm,test_sig,'test_ref','','','')

ref_sig = repmat(ref_temp,1,20);
ref_sig = ref_sig/max(abs(ref_sig));


% == Segmentation using ECGPUWAVE
[~,config] = wfdbloadlib;
ecgpuwave(recName,'edr',[],[],'qrs'); % important to specify the QRS because it seems that ecgpuwave is crashing sometimes otherwise
[allanns,type] = rdann(recName,'edr'); % read the annotation from ecgpuwave
disp('Ah!')

