function sqi = mcsqi2(raw,fecg,mecg,fs)
%mcSQI Spectral Coherence metric
% 
% Returns the coherence SQI.
% 
% 
% Input:
%   raw:         Raw data segment
%   fecg:        Residual data segment
%   mecg:        Chest lead segment
% 
% Output:
%   sqi:            resulting sSQI for segment
% 
% Fetal Extraction Toolbox, version 1.0, February 2014
% Released under the GNU General Public License
%
% Copyright (C) 2014 Fernando Andreotti
% Dresden University of Technology, Institute of Biomedical Engineering
% fernando.andreotti@mailbox.tu-dresden.de
%
% Last updated : 09-03-2014
%
%
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
    

%[Prm,f]=mscohere(raw,mecg,[],[],1024,fs);
[Prf,f]=mscohere(raw,fecg,[],[],1024,fs);
%[Pmf,f]=mscohere(mecg,fecg,[],[],1024,fs);


sqi1 = mean(Prf(f>=60&f<=100)); % average coherence for 60-100Hz (no high frequent noise added)
%sqi2 = 1-mean(Pmf(f<=100));     % 1- average cross-spectrum MECG-FECG for f<100Hz

sqi = sqi1;


end



