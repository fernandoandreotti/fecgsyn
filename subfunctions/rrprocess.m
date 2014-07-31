function rr = rrprocess(n,hrmean,lfhfr,hrstd,sfrr,flo,fhi,flostd,fhistd)
% this function simulates realistic hr time series
%
% inputs
%   N:      approximate number of heart beats [int]
%   hr:     mean heart rate [bpm]
%   lfhfr:  low to high freauency ratio of the two Gaussians
%           hf->simulates respiratory sinus arrythmia
%           lf->simulates Mayer
%   hrstd:  standard deviation of heart rate [bpm]
%   flo:    center freqency of low frequency Gaussian (Mayer) [Hz]
%   flhi:   center freqency of high frequency Gaussian (RSA) [Hz]
%
% output
%   rr: rr interval
%
% Last updated by Joachim Behar 23-12-2013 
%
% reference
% [1] McSharry, Patrick E., et al. "A dynamical model for generating 
% synthetic electrocardiogram signals." Biomedical Engineering, 
% IEEE Transactions on 50.3 (2003): 289-294.
%
%
% Copyright (c) 2003 by Patrick McSharry & Gari Clifford, All Rights Reserved  
% See IEEE Transactions On Biomedical Engineering, 50(3), 289-294, March 2003. and
% Advanced tools for ECG Signal Processing, Chapter 4 (Models) - for updates 
% Contact P. McSharry (patrick AT mcsharry DOT net) or 
% G.D. Clifford (gari AT alum DOT mit DOT edu)
%
%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 2 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program; if not, write to the Free Software
%    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
% 
% ecgsyn.m and its dependents/variants are freely availble from Physionet - 
% http://www.physionet.org/ - please report any bugs to the authors above.

% set parameter default values
if isempty(n)
   n = 512;
end
if isempty(hrmean)
   hrmean = 60;
end
if isempty(lfhfr)
   lfhfr = 0.5;
end
if isempty(hrstd)
   hrstd = 1;
end
if isempty(sfrr)
   sfrr = 1;
end
if isempty(flo)
   flo = 0.1;
end
if isempty(fhi)
   fhi = 0.25;
end
if isempty(flostd)
  flostd = 0.01;
end
if isempty(fhistd)
  fhistd = 0.01;
end


w1 = 2*pi*flo;
w2 = 2*pi*fhi;
c1 = 2*pi*flostd;
c2 = 2*pi*fhistd;
sig2 = 1;
sig1 = lfhfr;
rrmean = 60/hrmean;
rrstd = 60*hrstd/(hrmean*hrmean);

df = sfrr/n;
w = (0:n-1)'*2*pi*df;
dw1 = w-w1;
dw2 = w-w2;

Hw1 = sig1*exp(-0.5*(dw1/c1).^2)/sqrt(2*pi*c1^2);
Hw2 = sig2*exp(-0.5*(dw2/c2).^2)/sqrt(2*pi*c2^2);
Hw = Hw1 + Hw2;
Hw0 = [Hw(1:floor(n/2)); Hw(floor(n/2):-1:1)];
Sw = (sfrr/2)*sqrt(Hw0);

ph0 = 2*pi*rand(floor(n/2)-1,1);
ph = [ 0; ph0; 0; -flipud(ph0) ]; 
SwC = Sw .* exp(1j*ph);
x = (1/n)*real(ifft(SwC));

xstd = std(x);
ratio = rrstd/xstd;
rr = rrmean + x*ratio;

