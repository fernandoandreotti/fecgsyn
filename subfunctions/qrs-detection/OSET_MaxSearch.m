% peaks = OSET_MaxSearch(x,f,flag),
% R-peak detector based on max search
%
% inputs:
% x_ini: vector of input data
% ff: approximate ECG beat-rate in Hertz, normalized by the sampling frequency
% flag: search for positive (flag=1) or negative (flag=0) peaks. By default
% the maximum absolute value of the signal, determines the peak sign.
%
% output:
% peaks: vector of R-peak impulse train
%
% Notes:
% - The R-peaks are found from a peak search in windows of length N; where 
% N corresponds to the R-peak period calculated from the given f. R-peaks 
% with periods smaller than N/2 or greater than N are not detected.
% - The signal baseline wander is recommended to be removed before the
% R-peak detection
%
%
% Open Source ECG Toolbox, version 1.0, November 2006 
% Released under the GNU General Public License
% Copyright (C) 2006  Reza Sameni
% Sharif University of Technology, Tehran, Iran -- GIPSA-Lab, INPG, Grenoble, France
% reza.sameni@gmail.com
%
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
%
% Modifications by Fernando Andreotti (October 2013)
%         - Made selection of maxima/minima automatic and more robust
function peaks = OSET_MaxSearch(x_ini,ff,varargin)

if size(x_ini,1) > size(x_ini,2)
    x_ini = x_ini';
end
loopch = size(x_ini,1);  % attributes number of channels N to loopch
peaksmat = zeros(size(x_ini));

if(nargin==3)
    flag = varargin{1};
else
    % By choosing median of 5% of the data length makes flag choice more robust.
    % modification by Andreotti
    for k = 1: loopch
        x = decimate(x_ini(k,:),5);    % just to reduce calculus
        y = sort(x);    
        flag(k) = mean(abs(y(round(end-0.1*length(y)):end)))>mean(abs(y(1:round(0.1*length(y)))));
    end
end
% loops through every channel
for i=1:loopch                
    x = x_ini(i,:);
    N = length(x);
    peaks = zeros(1,N);
    
    th = .5;
    rng = floor(th/ff);
       
    if(flag)
        for j = 1:N,
            %         index = max(j-rng,1):min(j+rng,N);
            if(j>rng && j<N-rng)
                index = j-rng:j+rng;
            elseif(j>rng)
                index = N-2*rng:N;
            else
                index = 1:2*rng;
            end
            
            if(max(x(index))==x(j))
                peaks(j) = 1;
            end
        end
    else
        for j = 1:N,
            %         index = max(j-rng,1):min(j+rng,N);
            if(j>rng && j<N-rng)
                index = j-rng:j+rng;
            elseif(j>rng)
                index = N-2*rng:N;
            else
                index = 1:2*rng;
            end
            
            if(min(x(index))==x(j))
                peaks(j) = 1;
            end
        end
    end
    % remove fake peaks
    I = find(peaks);
    d = diff(I);
    % z = find(d<rng);
    peaks(I(d<rng))=0;
    peaks = find(peaks);
end
