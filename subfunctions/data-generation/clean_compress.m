function out=clean_compress(out)
% function out=clean_compress(out)
% this function eliminates some of the substructures from "out" and
% compresses the variables to int16 for saving disk space before saving
% 
% Input:
%    internal structure "out"
% 
%  Output:
%         same as input
%
% Examples:
% TODO
%
% See also:
% run_ecg_generator
% 
% 
% --
% fecgsyn toolbox, version 1.2, March 2017
% Released under the GNU General Public License
%
% Copyright (C) 2017  Joachim Behar & Fernando Andreotti
% Department of Engineering Science, University of Oxford
% joachim.behar@oxfordalumni.org, fernando.andreotti@eng.ox.ac.uk
%
% 
% For more information visit: http://www.fecgsyn.com
% 
% Referencing this work
%
% Behar, J., Andreotti, F., Zaunseder, S., Li, Q., Oster, J., & Clifford, G. D. (2014). An ECG Model for Simulating 
% Maternal-Foetal Activity Mixtures on Abdominal ECG Recordings. Physiol. Meas., 35(8), 1537–1550.
% 
%
% Last updated : 15-03-2017
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
%


    gain = 3000;
    out_tmp=rmfield(out,{'f_model' 'm_model' 'vols' 'selvcgm' 'selvcgf'});
    out = struct();
    out.mecg = int16(round(3000*out_tmp.mecg));
    if ~isempty(out_tmp.fecg)
        for i = 1:length(out_tmp.fecg)
            out.fecg{i} = int16(round(3000*out_tmp.fecg{i}));
        end
    else
        out.fecg = {};
    end
    if ~isempty(out_tmp.noise)
        for i = 1:length(out_tmp.noise)
            out.noise{i} = int16(round(gain*out_tmp.noise{i}));
        end
    else
        out.noise = {};
    end
    out.mqrs = out_tmp.mqrs;
    out.fqrs = out_tmp.fqrs;
    out.param = out_tmp.param;
end