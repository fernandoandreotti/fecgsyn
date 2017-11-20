function [out_comps,qrsmethod,W] = FECGSYN_bss_extraction(data,method,fs,blen,defl)
% Uses Blind Source Separation Methods for FECG extraction given a
% reference QRS. The component which most agrees with the reference, in
% terms of F1-measure, is picked as best channel.
%
% Available methods:
% Independent Component Analysis (ICA)
% Principal Component Analysis (PCA)
%
% Input
% data:      Matrix containing signals to serve as input for bss technique
% method:    String containing method name i.e. 'ICA' or 'PCA'
% fs:        Sampling frequency [Hz]
% blen:      Divide signal into segment of blen length
% defl:      Boolean, if 1 will apply PCA deflation to data
% refqrs:    Array containing reference QRS detections for F1 measure
% (optional)
% blen:      Iterates method every blen (in seconds)
% filename:  Saves output of ica into filename
%
% Output
% out_comps:    selected best channels on every block
% qrsmethod:    BSS techniques chosen QRS detections   
% W:            used mixing matrices
%
%
% Examples:
% TODO
%
% See also:
% FECGSYN_ts_extraction
% FECGSYN_kf_extraction
% FECGSYN_adaptfilt_extraction
% FEGSYN_main_extract
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

%% == Input test
if size(data,1) > size(data,2)
    data = data';
end

method = upper(method);
Bold = [];	% old mix matrix

%% Main
% Loop BSS technique on every block
qrsmethod = [];
blen = blen * fs;
ssamp = 1; % starting sample to filter (offset)
endsamp = ssamp + blen - 1;      % ending sample to filter
loop = 1; % allows iterations

if size(data,2)<endsamp
    blen = size(data,2);
    endsamp = size(data,2);
end


W = cell(5,1);
count = 0;
out_comps = {};  % allocating

while (loop)  % will quit as soon as complete signal is filtered
     if (size(data,2) - ssamp) < 1.5*blen     % if there is less than 1.5x blen
        endsamp = size(data,2);              % interval, it should filter until end
        loop = 0;                            % and be the last loop iteration
     end
    count = count +1;
    samp2filt = ssamp:endsamp;              % creating a list with samples to filter
    if defl 
    % = use PCA in keeping channels with 95% of eigenspectrum
     % = normalising data (PCA is sensible to scaling)
     tmpdata = data(:,samp2filt);
     tmpdata = bsxfun(@minus,tmpdata,mean(tmpdata,2)); % remove mean (JADE is sensible)
     Nmat = diag(std(tmpdata,0,2));   % normalization matrix
     tmpdata = bsxfun(@rdivide,tmpdata,std(tmpdata,0,2)); % divide by standard deviation
     [coeff, score, latent] = pca(tmpdata');
     perc = cumsum(latent)./sum(latent);
     Ncomp = find(perc>=0.999,1,'first');   % keeping 99.9% data variance
     tmpdata = score(:,1:Ncomp)*coeff(:,1:Ncomp)';
     tmpdata = tmpdata';
     
    else
        tmpdata = data(:,samp2filt);
        tmpdata = bsxfun(@minus,tmpdata,mean(tmpdata,2)); % remove mean (JADE is sensible)
        Nmat = diag(std(tmpdata,0,2));   % normalization matrix
        tmpdata = bsxfun(@rdivide,tmpdata,std(tmpdata,0,2)); % divide by standard deviation
    end
          
    % this is because FastICA is not deterministic so make sure to use the same random seed at 
    % each run
    stream = RandStream.getGlobalStream; 
    reset(stream);
    switch method
        case 'FASTICA_DEF'
            % FastICA with deflation appraoch
            [~,~,Bnew] = fastica(tmpdata,'g','tanh','verbose','on','maxNumIterations',size(tmpdata,1)*1000,'approach','defl');
            disp(['FASTICA_DEF output size:' num2str(size(Bnew,1)) 'x' num2str(size(Bnew,2))])
        case 'FASTICA_SYM'
            % FastICA with symmetric method
            [~,~,Bnew] = fastica(tmpdata,'g','tanh','verbose','on','maxNumIterations',size(tmpdata,1)*1000,'approach','symm');    
            disp(['FASTICA_SYM output size:' num2str(size(Bnew,1)) 'x' num2str(size(Bnew,2))])
        case 'JADEICA'
            % JADEICA (no restriction on number of sources)
            Bnew = jadeR(tmpdata,Ncomp); 
        case 'PCA'
            [Bnew,~] = pca(tmpdata');
            Bnew = Bnew';            
        otherwise
            error('bss_extraction: Method not implemented')
    end
    if isempty(Bold); Bold = Bnew; end; % first loop
    outnew = Bold*Nmat*tmpdata;
    outnew = diag(1./max(outnew,[],2))*outnew; % may not have the same size as "data"
    W{count} = Bold;   % saving previous mixing matrices
    Bold = Bnew;   
    
    % Adding to previous blocks, differences in the number of components
    % are corrected with zeros
    out_comps{count} = outnew;
    

    %Augment offsets
    ssamp = endsamp+1;
    endsamp = endsamp+blen;
end

pad = max(cellfun(@(x) size(x,1),out_comps));
for k = 1:length(out_comps)
    out_comps{k} = [out_comps{k}; zeros(pad-size(out_comps{k},1),length(out_comps{k}))];
end
out_comps = cell2mat(out_comps);



