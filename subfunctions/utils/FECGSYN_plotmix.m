function FECGSYN_plotmix(out)
% Function to plot generated mixtures
% 
% This functin plots Nplots channels from fecgsyn's internal struct "out".
% It is useful in checking the results of data generation during debug.
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
% Referencing this work
%
%   Behar Joachim, Andreotti Fernando, Zaunseder Sebastian, Li Qiao, Oster Julien, Clifford Gari D. 
%   An ECG simulator for generating maternal-foetal activity mixtures on abdominal ECG recordings. 
%   Physiological Measurement.35 1537-1550. 2014.
% 
% 
% --
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

FONT_SIZE = 15;
LINE_WIDTH = 1.2;
MARKER_SIZE = 5;
Nplots = 5;      % number of plots

if ~exist('out.mixture','var')
    if ~isempty(out.noise)
        out.mixture = double(out.mecg) + sum(cat(3,out.fecg{:}),3) ...
            + sum(cat(3,out.noise{:}),3);     % re-creating abdominal mixture
    else
        out.mixture = double(out.mecg) + sum(cat(3,out.fecg{:}),3);
    end
    
end
Nchan = size(out.mixture,1);
chan = randsample(1:Nchan,Nplots);
figure(1)
clf

for i = 1:Nplots
    h(i) = subplot(Nplots,1,i);   
    plot(out.mixture(chan(i),:),'Color',[.6 .6 .6],'LineWidth',LINE_WIDTH)
    hold on
    for j = 1:length(out.fqrs)
        valf = 1.2*median(out.mixture(chan(j),out.fqrs{j}))*ones(1,length(out.fqrs{j}));
        plot(out.fqrs{j},valf,'xr','MarkerSize',MARKER_SIZE,'MarkerFaceColor','r')
    end
    valm = 1.2*median(out.mixture(chan(i),out.mqrs))*ones(1,length(out.mqrs));
    plot(out.mqrs,valm,'dk','MarkerSize',MARKER_SIZE,'MarkerFaceColor','k')
    ylabel(['ch' num2str(chan(i))],'FontSize',FONT_SIZE)
    hold off
end
linkaxes(h)
