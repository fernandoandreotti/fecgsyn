function [outsig,qrsmethod] = FECGSYN_bss_extraction(data,method,fs,refqrs,varargin)
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
% refqrs:    Array containing reference QRS detections for F1 measure
% (optional)
% blen:      Iterates method every blen (in seconds)
% filename:  Saves output of ica into filename
%
% Output
% outsig:    selected best channels on every block
% qrsmethod: BSS techniques chosen QRS detections   
%
%
% NI-FECG simulator toolbox, version 1.0, February 2014
% Released under the GNU General Public License
%
% Copyright (C) 2014  Joachim Behar & Fernando Andreotti
% Oxford university, Intelligent Patient Monitoring Group - Oxford 2014
% joachim.behar@eng.ox.ac.uk, fernando.andreotti@mailbox.tu-dresden.de
%
% Last updated : 22-07-2014
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

%% == Input test
if size(data,1) > size(data,2)
    data = data';
end

switch length(varargin)
    case 0
        blen = size(data,2);
    case 1
        blen = varargin{1};
    case 2
        blen = varargin{1};
        outfilename = varargin{2};
    otherwise
        error('ica_extraction: too many inputs given to function')
end
method = upper(method);

% Parameters for QRS detector
INTERV = round(0.05*fs);     % BxB acceptance interval
REFRAC = .15; % detector refractory period (s)
TH = 0.3; 	% detector threshold
Bold = [];	% old mix matrix

%% Main
% Loop BSS technique on every block
qrsmethod = [];
blen = blen * fs;
ssamp = 1; % starting sample to filter (offset)
endsamp = ssamp + blen - 1;      % ending sample to filter
loop = 1; % allows iterations

if size(data,2)<endsamp
    blen = Data.length;
    endsamp = size(data,2);
end


while (loop)  % will quit as soon as complete signal is filtered
    if (size(data,2) - ssamp) < 1.5*blen     % if there is less than 1.5x blen
        endsamp = size(data,2);              % interval, it should filter until end
        loop = 0;                            % and be the last loop iteration
    end
    
    samp2filt = ssamp:endsamp;              % creating a list with samples to filter
    idx = refqrs>=ssamp & refqrs <= endsamp;
    refint = refqrs(idx) - ssamp; % within interval
    % this is because FastICA is not deterministic so make sure to use the same random seed at 
    % each run
    stream = RandStream.getGlobalStream; 
    reset(stream);
    switch method
        case 'FASTICA_DEF'
            % FastICA with deflation appraoch
            [~,~,Bnew] = fastica(data(:,samp2filt),'g','tanh','verbose','on','maxNumIterations',size(data,1)*1000,'approach','defl');
            disp(['FASTICA_DEF output size:' num2str(size(Bnew,1)) 'x' num2str(size(Bnew,2))])
        case 'FASTICA_SYM'
            % FastICA with symmetric method
            [~,~,Bnew] = fastica(data(:,samp2filt),'g','tanh','verbose','on','maxNumIterations',size(data,1)*1000,'approach','symm');    
            disp(['FASTICA_SYM output size:' num2str(size(Bnew,1)) 'x' num2str(size(Bnew,2))])
        case 'JADEICA'
            % JADEICA (no restriction on number of sources)
            Bnew = jadeR(data(:,samp2filt)); 
        case 'PCA'
            [Bnew,~] = princomp(data(:,samp2filt)');
            Bnew = Bnew';
        otherwise
            error('bss_extraction: Method not implemented')
    end
    if isempty(Bold); Bold = Bnew; end;
    outdata = Bold*data(:,samp2filt);
    outdata = diag(1./max(outdata,[],2))*outdata; % may not have the same size as "data"
    Bold = Bnew;   

    % = QRS detect each component and take F1, RMS measure
    qrsdet = cell(1,size(outdata,1));
    F1 = zeros(1,size(outdata,1));
    RMS = zeros(1,size(outdata,1));
    for ch = 1:size(outdata,1)
        qrsdet{ch} = qrs_detect(outdata(ch,:),TH,REFRAC,fs);
        if ~isempty(qrsdet{ch})
            [F1(ch),RMS(ch)] = Bxb_compare(refint,qrsdet{ch},INTERV);
        else
            F1(ch) = 0;
            RMS(ch) = NaN;
        end
    end
    
    % = decide for one component
    %Saving segment
    [~,maxch] = max(F1);
    qrsmethod = [qrsmethod (qrsdet{maxch}+ssamp-1)];
    outsig(samp2filt) = outdata(maxch,:);
    
    if length(varargin)==2
        if ssamp > 1
            outdatanew = outdata;
            maxchnew = maxch;
            load([outfilename '_' method])           
            % treating cases when BSS technique outputs less channels than before
            if size(outdatanew,1)<size(outdata,1)
                outdatanew = [outdatanew;zeros(size(outdata,1)-size(outdatanew,1),length(outdatanew))];
                outdata = [outdata outdatanew];
            elseif size(outdatanew,1)>size(outdata,1) % case first run outputed less components
                outdata = [outdata;zeros(size(outdatanew,1)-size(outdata,1),length(outdata))];
                outdata = [outdata outdatanew];
            else % trivial case
                outdata = [outdata outdatanew];
            end
            maxch = [maxch maxchnew];
        end
        save([outfilename '_' method],'outdata','maxch');
    end
    
    %Augment offsets
    ssamp = endsamp+1;
    endsamp = endsamp+blen;
    

end


