function[qrs_pos,varargout] = pantompkins_qrs(data,fs,varargin)
%function [qrs_pos] = pantompkins_qrs(data,fs,varargin)
%
% QRS Detection with Pan-Tompkins algorithm. 
% Pan et al. 1985: A Real-Time QRS Detection Algorithm.
% IEEE Trans Bio Eng. Vol. 32, No. 3, S.230-236
%
% This is an offline implementation of the detection, which nevertheless
% is realizing a adaptive detection (from beat to beat) like the proposed
% online implementation. Compared to the original proposed filtering, all 
% filtering is replaced by phase zero filtering.
%
% Impementation includes a powerline-bandstop-filter at 50 Hz, change this
% characteristic if necessary inside included function pantompkinsFilter().
% Furthermore, the widest possible QRS-Complex is initialized as having a
% duration of 150 ms. This could be adjusted inside included function
% pantompkinsFilter() as well. Another parameters are the refractory period
% and the t-wave exclusion period which are setted up inside the main
% code-function as refrac = 200 ms and refracT = 360 ms.
% Furthermore, the threshold-usage (as originally proposed, both thresholds
% of filtered and integrated data have to be exceeded for valid qrs) can be
% updated in terms of allowing one of both thresholds to deceed level, if an
% associated region is detected (for example a short threshold-deceeding,
% which nevertheless belongs to the same QRS-complex). This option is
% activated by regionFLAG = 1. 
%
%
% input: - data: input-data (single channel vector)
%        - fs: data sampling-frequency
%        - varargin
%        - verbose:   enable/disables function verbose [bool]
%        - str_ident: dataset-identification-string used for error-messages
% output:
%        - qrs_pos: indexes of detected peaks (in samples)
%        - varargout:
%        - {1}: filtered data (after powerline and bandpass filtering)
%        - {2}: integrated data (after complete preprocessing)
%        - {3}: threshold-array for filtered data (threshold for each QRS-complex)
%        - {4}: threshold-array for integrated data (threshold for each QRS-complex)
%
% Released under the GNU General Public License
%
% Copyright (C) 2014  Daniel Wedekind
% Technical University Dresden
% Faculty of Electrical and Computer Engineering
% Institute of Biomedical Engineering
% Cardiovascular Signal Processing Group
% Dresden, Germany, 2014
%
% daniel.wedekind@mailbox.tu-dresden.de
%
% Vers. 1.0:
% Last updated : 29-05-2014
% 
% Change log:
%  29-05-2014: Andreotti   added refractory times and qrs width as input
%                          parameter
%
% To be included in future (possible improvements):
% - QRS-replacement after RR-analysis (replacement functionality is already
% included in the code, it just has to be enabled each time using the
% variable replaceQRS)
% - The T-Wave critereon (one half of the maximal slope of the preceeding
% QRS complex might perform to strict. this parameter for example could be
% updated (softened) if RR-irregularity has been detected, otherwise this
% could be combined with the above mentioned QRS-replacement which is
% already done but with weak QRS replacement power (a replacement option
% is only accessed in case of regular RR-criterion (92-116% RR-variance)
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


% Input test
optargs = {false,'DEFAULT',round(0.2*fs),round(.36*fs),round(.15*fs)};     % default values
newVals = cellfun(@(x) ~isempty(x), varargin);
optargs(newVals) = varargin(newVals);
[verbose,str_ident,refrac,refracT,QRSmaxwidth] = optargs{:};

if nargin > 7
    error('Wrong number of input arguments specified!');
end
if ~islogical(verbose)
    error('Verbose option must be boolean.')
end
if (size(data,1) > size(data,2))
    data = data';
elseif (size(data,1) == size(data,2))
    error('check input data dimension. input data vector is of size 1x1');
elseif (~((size(data,1) == 1)||(size(data,2) == 1)))
    error('check input data dimension. input data vector seems to contain multiple channels');
end

if (size(data,2) < (3*fs))
    error('check input data dimension. input data vector containing less then 3s of data');
end

qrs_pos = [];
thresholdI1m = [];
thresholdF1m = [];
%% allow one of both thresholds I or F to be deceeded if associated region is detected   
regionFLAG = 1;

%% general signal preparation
% baseline-wander
[BWb,BWa] = butter(5,[1.0].*2/fs,'high');
data = filtfilt(BWb,BWa,data);
% normalization
data = data./std(data);

%% treat case of negative QRS sign
checkSign = sort(data,'descend');
signFLAG = mean(abs(checkSign(length(data)-round(0.1*length(checkSign)):length(data))))>mean(abs(checkSign(1:round(0.1*length(checkSign)))));
if (signFLAG) && verbose        
%      disp('input channel inverted due to possible negative QRS manifestation');
     data = -data;
end

%% applying pan tompkins filtering
[filt_dat,diff_dat,mov_dat] = pantompkinsFilter(data, fs,QRSmaxwidth);

%% learning phase 1
% initializing thresholds based on first 2s of signal 
peakI = max(mov_dat(1,1:(2*fs)));
thresholdI1 = 0.125*peakI;
%thresholdI2 = 0.5*thresholdI1; % is automatically updated if necessary  
peakF = max(filt_dat(1,1:(2*fs)));
thresholdF1 = 0.125*peakF;
%thresholdF2 = 0.5*thresholdF1; % is automatically updated if necessary

%% learning phase 2 (this isn't explicitly specified in the originally proposed paper)
learned = 0;
RR_init = 0;
count_thresAbove = 0;
i = 1;

%including possibility to reduce threshold to one half for five times during learning phase 2
%determine threshold-exceeding signal excerpts
while ((learned == 0)&&(i <= 6))
flagS = 0;
flagE = 0;
    
aboveI1 = (mov_dat(1,:)>thresholdI1);
aboveF1 = (filt_dat(1,:)>thresholdF1);
%exceeding both thresholds
aboveCommon = aboveI1.*aboveF1;
negativeCommon = ~(aboveCommon);
%assign regions
aboveCommonRegion = aboveCommon;
found_region = 0;
start_index = 0;
for index = 1:length(aboveCommonRegion)
    %starting point of region
    if(aboveCommonRegion(1,index) == 1)
        found_region = 1;
        start_index = index;
    end
    %assign region
    if(found_region)
        aboveCommonRegion(1,index) = aboveI1(1,index)|aboveF1(1,index);
    %detect end of region
    regI1 = diff(aboveI1(1,start_index:index));
    regF1 = diff(aboveF1(1,start_index:index));
    if(~((sum(abs(regI1))==0)||(sum(abs(regF1))==0)))
        found_region = 0;
    end
    if((sum(abs(regI1))>0)&&(sum(abs(regF1))>0))
        aboveCommonRegion(1,start_index+1:index) = 0;
    end
    end
end

%count exceeding regions
aboveBackup = aboveCommon;
if (regionFLAG)
    aboveCommon = aboveCommonRegion;
end

count_thresAbove = sum(abs(diff(aboveCommon)));
if (aboveCommon(1,1) == 1)
    count_thresAbove = count_thresAbove + 1;
    flagS = 1;
end
if (aboveCommon(1,length(aboveCommon)) == 1)
    count_thresAbove = count_thresAbove + 1;
    flagE = 1;
end
count_thresAbove = count_thresAbove/2;
%threshold re-adjustment for possible learning phase recall (learned = 0)
thresholdI1 = thresholdI1/2;
thresholdF1 = thresholdF1/2;
i = i + 1;
if (i >= 3)&& verbose
    %disp(['QRS-detector notification (Pan&Tompkins):']);
    disp(['threshold reduced during learning phase of dataset' str_ident ]);
end

% generation of corresponding index-pairs for threshold-exceeding regions
start_region_index = (find(diff(aboveCommon) == 1))+1;
end_region_index = find(diff(aboveCommon) == -1);
if (flagS)
    start_region_index = [1,start_region_index];
end
if (flagE)
    end_region_index = [end_region_index, length(aboveCommon)];
end

j = 1;
k = 2;

%determine first two reliable initial qrs-complexes
while ((learned == 0)&&(j < count_thresAbove)&&(k <= count_thresAbove)&&(count_thresAbove >= 2))
    
    %pick out current threshold-exceeding signal regions
    first_region = data(start_region_index(j):end_region_index(j));
    second_region = data(start_region_index(k):end_region_index(k));
    
    %calculating distance of peaks using input signal 
    [cmax1,imax1] = max(first_region);
    imax1 = imax1+start_region_index(j)-1;
    [cmax2,imax2] = max(second_region); 
    imax2 = imax2+start_region_index(k)-1;
    distance = imax2-imax1;
    
    %check refractory/t-wave criterion and manipulate indices
    %inside refractory period
    if (distance <= refrac)
        k = k+1;
    %between refractory period and t-wave-citerion     
    elseif ((distance <= refracT)&&(distance > refrac))
        slope_ratio = max(diff_dat(start_region_index(j):end_region_index(j)))/max(diff_dat(start_region_index(k):end_region_index(k))); 
        if (slope_ratio > 2)
            %probably t-wave
            k = k+1;
        elseif (slope_ratio < 0.5)
            %j = j+1; %preceeding QRS-detection is identified as t-wave 
            %k = k+1;
            learned = 1;       
        else
            learned = 1;
        end
        
    else
        learned = 1;
    end
           
end

end

if (learned ~= 0)
    qrs_pos(1,1) = 0;
    qrs_pos(1,2) = imax1;
    qrs_pos(1,3) = imax2;
    %re-initialize thresholds for qrs-detection regarding first two detections
    RR_init = imax2-imax1;
    if verbose;disp(['RR_init: ' num2str(RR_init/fs) ' s']);end
    
    posMask = aboveBackup(start_region_index(j):end_region_index(k));
    negMask = negativeCommon(start_region_index(j):end_region_index(k));
    
    peakIs = max(posMask.*mov_dat(start_region_index(j):end_region_index(k)));
    peakIn = max(negMask.*mov_dat(start_region_index(j):end_region_index(k)));
    
    peakFs = max(posMask.*filt_dat(start_region_index(j):end_region_index(k)));
    peakFn = max(negMask.*filt_dat(start_region_index(j):end_region_index(k)));
   
    SPKI = peakIs;
    NPKI = peakIn;
    SPKF = peakFs;
    NPKF = peakFn;
    
    thresholdI1m(1,1) = thresholdI1*2;
    thresholdF1m(1,1) = thresholdF1*2;
    
    % [cmax, starting_detection] = max(posMask.*data(start_region_index(j):end_region_index(k)));
    % qrs_pos(1,1) = starting_detection+start_region_index(j)-1;
    qrs_pos(1,1) = imax1;
end

if (learned == 0)
    warning(['QRS-detection learning phase of dataset ' str_ident ' was not successful.']);
    return
end



%% detection phase
%initialization of detection phase
qrs_pos = qrs_pos(1,1);

RR_AV1vec = NaN(1,8);
RR_AV1vec(1,1) = RR_init;
RR1 = 0;

RR_AV2vec = NaN(1,8);
RR_AV2vec(1,1) = RR_init;
RR2 = 0;

finished = 0;
thres_add = 1;
end_add = 0;
segment_OK = 0;
initial_TA = 1;
proc_per_QRS = 1;
actual_index = qrs_pos(1,1);

while ((finished == 0)&&(ceil((1.66*RR1)+actual_index)<=length(data)))

    if(segment_OK)
        end_add = 0;
        thres_add = 1;
        initial_TA = 1;
    end
    
    RR1 = mean(RR_AV1vec(find(~isnan(RR_AV1vec))));
    RR2 = mean(RR_AV2vec(find(~isnan(RR_AV2vec))));
    
    %regulary heart rate check
    if (RR1 == RR2)
        factor = 1;
    else
        factor = 0.5;
    end    
    
    %updating thresholds 
    thresholdI1 = factor*(NPKI + (0.25*(SPKI-NPKI)))*thres_add;
    %thresholdI2 = 0.5*thresholdI1; %is automatically updated if necessary

    thresholdF1 = factor*(NPKF + (0.25*(SPKF-NPKF)))*thres_add;
    %thresholdF2 = 0.5*thresholdF1; %is automatically updated if necessary
    
    %begin new search from actual detection
    %building up data and apply threshold
    if ((actual_index+(1.66*RR1)+end_add) <= length(data))
    end_segment = ceil((actual_index+(1.66*RR1)))+end_add;
    actual_mov_dat = mov_dat(actual_index:end_segment);
    actual_filt_dat = filt_dat(actual_index:end_segment);
    actual_diff_dat = diff_dat(actual_index:end_segment);
    
    else
        break;
    end
    
    i = 1;
    segment_OK = 0;
    
    
    %including possibility to reduce threshold to one half for one time
    %during detection phase of one single segment 
    %determine threshold-exceeding signal excerpts
    while((segment_OK == 0)&&(i <= 2))
        
        flagS = 0;
        flagE = 0;
    
        aboveI1 = (actual_mov_dat(1,:)>thresholdI1);
        aboveF1 = (actual_filt_dat(1,:)>thresholdF1);
        %special case of every feature signal value is below threshold after threshold-halving
        if (min(actual_mov_dat(1,:)) > thresholdI1)
            
            actual_range = max(actual_mov_dat(1,:)) - min(actual_mov_dat(1,:));
            
            new_thresM = exp(-(proc_per_QRS/3));
            new_thresM = new_thresM*actual_range;
            new_thresM = min(actual_mov_dat(1,:))+new_thresM;
            
            aboveI1 = (actual_mov_dat(1,:) > new_thresM);
            
        end
        if (min(actual_filt_dat(1,:)) > thresholdF1)
            
            actual_range = max(actual_filt_dat(1,:)) - min(actual_filt_dat(1,:));
            
            new_thresF = exp(-(proc_pre_QRS/3));
            new_thresF = new_thresF*actual_range;
            new_thresF = min(actual_filt_dat(1,:))+new_thresF;
            
            aboveF1 = (actual_filt_dat(1,:) > new_thresF);
                       
        end
        
        
        %exceeding both thresholds
        aboveCommon = aboveI1.*aboveF1; 
        negativeCommon = ~(aboveCommon);
        %assign regions
        aboveCommonRegion = aboveCommon;
        found_region = 0;
        start_index = 0;
        for index = 1:length(aboveCommonRegion)
        %starting point of region
        if(aboveCommonRegion(1,index) == 1)
            found_region = 1;
            start_index = index;
        end
        %assign region
        if(found_region)
            aboveCommonRegion(1,index) = aboveI1(1,index)|aboveF1(1,index);
        %detect end of region
        regI1 = diff(aboveI1(1,start_index:index));
        regF1 = diff(aboveF1(1,start_index:index));
        if(~((sum(abs(regI1))==0)||(sum(abs(regF1))==0)))
            found_region = 0;
        end
        if((sum(abs(regI1))>0)&&(sum(abs(regF1))>0))
            aboveCommonRegion(1,start_index+1:index) = 0;
        end
        end
        end
        
        %treating case of beeing unable to detect most recent detection again after updating threshold
        if (aboveCommon(1,1) == 0)
            aboveCommon(1,1) = 1;
            aboveCommonRegion(1,1) = 1;
        end 
        
        
        %count exceeding regions
        aboveBackup = aboveCommon;
        if (regionFLAG)
            aboveCommon = aboveCommonRegion;
        end
                
        count_thresAbove = sum(abs(diff(aboveCommon)));
        if (aboveCommon(1,1) == 1)
            count_thresAbove = count_thresAbove + 1;
            flagS = 1;
        end
        if (aboveCommon(1,length(aboveCommon)) == 1)
            count_thresAbove = count_thresAbove + 1;
            flagE = 1;
        end
        count_thresAbove = count_thresAbove/2;
        %threshold re-adjustment for possible searchback (segment_OK = 0)
        thresholdI1 = thresholdI1/2;
        thresholdF1 = thresholdF1/2;
        i = i + 1;
        if (i >= 3)
            proc_per_QRS = proc_per_QRS+1;
            if (verbose)
            disp(['threshold reduced during detection phase of dataset ' str_ident ' - segment ' num2str(actual_index) '-' num2str(end_segment)]);
            end
        end

        % generation of corresponding index-pairs for threshold-exceeding regions
        start_region_index = (find(diff(aboveCommon) == 1))+1;
        end_region_index = find(diff(aboveCommon) == -1);
        if (flagS)
            start_region_index = [1,start_region_index];
        end
        if (flagE)
            end_region_index = [end_region_index, length(aboveCommon)];
        end

        j = 1;
        k = 2;
        replaceQRS = 0;
        %determine the next reasonable QRS complex
        while ((segment_OK == 0)&&(k <= count_thresAbove)&&(count_thresAbove >= 2))
    
        %pick out current threshold-exceeding signal regions
        first_region = data((start_region_index(j)+actual_index-1):(end_region_index(j)+actual_index-1));
        second_region = data((start_region_index(k)+actual_index-1):(end_region_index(k)+actual_index-1));
    
        %calculating distance of peaks using input signal 
        [cmax1,imax1] = max(first_region);
        imax1 = imax1+start_region_index(j)-1;
        [cmax2,imax2] = max(second_region); 
        imax2 = imax2+start_region_index(k)-1;
        distance = imax2-imax1;
    
        %check refractory/t-wave criterion and manipulate indices
        %inside refractory period
        if (distance <= refrac)
            
            aboveCommon(start_region_index(k):end_region_index(k)) = 0; %due to subsequent threshold adjustment 
            negativeCommon(start_region_index(k):end_region_index(k)) = 1;
            
            k = k+1;
            
        %between refractory period and t-wave-criterion     
        elseif ((distance <= refracT)&&(distance > refrac))
            slope_ratio = max(actual_diff_dat(start_region_index(j):end_region_index(j)))/max(actual_diff_dat(start_region_index(k):end_region_index(k))); 
            if (slope_ratio > 2)
                %probably t-wave
                
                aboveCommon(start_region_index(k):end_region_index(k)) = 0; %due to subsequent threshold adjustment
                negativeCommon(start_region_index(k):end_region_index(k)) = 1;
                
                k = k+1;
                
            elseif (slope_ratio < 0.5)
                %replaceQRS = 1; %preceeding QRS-detection is identified as t-wave and replaced
                
                segment_OK = 1;
                end_add = 0;
                thres_add = 1;
                proc_per_QRS = 1;
            else
                segment_OK = 1;
                end_add = 0;
                thres_add = 1;
                proc_per_QRS = 1;
            end
            
        else
            segment_OK = 1;
            end_add = 0;
            thres_add = 1;
            proc_per_QRS = 1;
            
            %check if preceeding detection could be a t-wave
            if (length(qrs_pos) >= 2)
                rr_new = (imax2+actual_index-1)-qrs_pos(1,length(qrs_pos));
                rr_old = qrs_pos(1,length(qrs_pos))-qrs_pos(1,length(qrs_pos)-1);
                comb_detection = rr_new+rr_old;
                if ((comb_detection >=  0.92*RR2)&&(comb_detection <= 1.16*RR2))
                    replaceQRS = 1;
                end
            end
            
        end
        
        end
    
    end
        
    if(segment_OK)
    %estimate detection inside input data
    detection = imax2+actual_index-1;    
    
    %sort in new detection
    if (replaceQRS == 0)
       qrs_pos(1,1+length(qrs_pos)) = detection;
       
       thresholdI1m(1,1+length(thresholdI1m)) = thresholdI1*2;
       thresholdF1m(1,1+length(thresholdF1m)) = thresholdF1*2;
    else
       qrs_pos(1,length(qrs_pos)) = detection;
       
       thresholdI1m(1,length(thresholdI1m)) = thresholdI1*2;
       thresholdF1m(1,length(thresholdF1m)) = thresholdF1*2;
    end
       
    
    %calculate and sort in new RR-interval 
    if (length(qrs_pos) >= 2)
    RR_new = qrs_pos(1,length(qrs_pos))-qrs_pos(1,length(qrs_pos)-1);
        
    if (replaceQRS == 0)&& verbose
        disp(['recent RR-interval: ' num2str(RR_new/fs) ' s']);
    elseif verbose
        disp(['recent RR-interval: ' num2str(RR_new/fs) ' s (replaced)']);
    end
    
    if (isempty(find(isnan(RR_AV1vec(1,:)))))
        if(replaceQRS == 0)
            RR_AV1vec = RR_AV1vec(1,2:length(RR_AV1vec));
            RR_AV1vec(1,1+length(RR_AV1vec)) = RR_new;
        else
            RR_AV1vec(1,length(RR_AV1vec)) = RR_new;
        end
    else
        index_insert = find(isnan(RR_AV1vec(1,:)));
        if(replaceQRS == 0)
            index_insert = index_insert(1,1);
            RR_AV1vec(1,index_insert) = RR_new;
        else
            index_insert = index_insert(1,1)-1;
            RR_AV1vec(1,index_insert) = RR_new;
        end
    end
    
    if ((RR_new >= (0.92*RR1))&&(RR_new <= (1.16*RR1)))
    
        if (isempty(find(isnan(RR_AV2vec(1,:)))))
            if(replaceQRS == 0)
                RR_AV2vec = RR_AV2vec(1,2:length(RR_AV2vec));
                RR_AV2vec(1,1+length(RR_AV2vec)) = RR_new;
            else
                RR_AV2vec(1,length(RR_AV2vec)) = RR_new;
            end
        else
            index_insert = find(isnan(RR_AV2vec(1,:)));
            if(replaceQRS == 0)
                index_insert = index_insert(1,1);
                RR_AV2vec(1,index_insert) = RR_new;
            else
                index_insert = index_insert(1,1)-1;
                RR_AV2vec(1,index_insert) = RR_new;
            end
        end
    elseif verbose
        disp('RR-regularity dumped - out of 92-116% RR-variance');
    end
    end    
    
    %update actual start/end index for next treated segment
    actual_index = detection;
    
    %update running estimates of peak-values regarding threshold usage
    if (replaceQRS == 0)
        posMask = aboveBackup(start_region_index(j):end_region_index(k));
    else
        posMask = aboveBackup(start_region_index(k):end_region_index(k));
    end
        negMask = negativeCommon(start_region_index(j):end_region_index(k));
        
    if (replaceQRS == 0)
        peakIs = max(posMask.*actual_mov_dat(start_region_index(j):end_region_index(k)));
        peakFs = max(posMask.*actual_filt_dat(start_region_index(j):end_region_index(k)));
    else
        peakIs = max(posMask.*actual_mov_dat(start_region_index(k):end_region_index(k)));
        peakFs = max(posMask.*actual_filt_dat(start_region_index(k):end_region_index(k)));
    end
    peakIn = max(negMask.*actual_mov_dat(start_region_index(j):end_region_index(k)));
    peakFn = max(negMask.*actual_filt_dat(start_region_index(j):end_region_index(k)));
        
   
   
    if (i <= 2)
        SPKI = (0.125*peakIs)+(0.875*SPKI);
        NPKI = (0.125*peakIn)+(0.875*NPKI);
        SPKF = (0.125*peakFs)+(0.875*SPKF);
        NPKF = (0.125*peakFn)+(0.875*NPKF);
    else %double speed threshold adjustment after searchback
        SPKI = (0.25*peakIs)+(0.75*SPKI);
        NPKI = (0.25*peakIn)+(0.75*NPKI);
        SPKF = (0.25*peakFs)+(0.75*SPKF);
        NPKF = (0.25*peakFn)+(0.75*NPKF);
    end
    
    end
    
    %handling case of no detection after applying both thresholds (this isn't explicitly specified in the originally proposed paper)
    
    if (segment_OK == 0)
        proc_per_QRS = proc_per_QRS + 1;
    if ((end_segment + end_add + fs) <= length(data))
        %adding 1s of remaining signal to the already investigated one
        end_add = end_add + fs;
        %reducing thresholds again
        thres_add = thres_add*0.5;
        if verbose;disp(['additive calculation loop due to missing QRS-detection after regular processing phase...']);end
        %compensation of threshold recalculation 
        if (initial_TA)
            thres_add = thres_add*0.5;
            initial_TA = 0;
        end
    else
        finished = 1;
    end
    end
end

varargout{1} = filt_dat;
varargout{2} = mov_dat;
varargout{3} = thresholdF1m;
varargout{4} = thresholdI1m;

return



function [filt_dat,diff_dat,mov_dat] = pantompkinsFilter(data, fs,QRSmaxwidth)
% Filtering of data matrix according to proposed data
% processing of QRS-detection algorithm:
% Pan et al. 1985: A Real-Time QRS Detection Algorithm.
% IEEE Trans Bio Eng. Vol. 32, No. 3, S.230-236
%
% Compared to the original proposed filtering, all filtering is replaced by 
% phase zero filtering.
%
% Impementation includes a powerline-bandstop-filter at 50 Hz, change this
% characteristic if necessary. Furthermore, the widest possible QRS-Complex 
% is initialized as having a duration of 150 ms. This could be adjusted as well. 
%
%
% input: - data: input-data of size MxN
%                M: amout of single input channels
%                N: single data channel length
%        - fs: data sampling-frequency
%        - varargin
%        NONE
% output:
%        - filt_dat: data after initial powerline and bandpass filtering
%        - diff_dat: filt_dat after smoothed differentation
%        - mov_dat:  diff_dat after moving window integration
%
% Released under the GNU General Public License
%
% Copyright (C) 2014  Daniel Wedekind
% Technichal University Dresden
% Faculty of Electrical and Computer Engineering
% Institute of Biomedical Engineering
% Cardiovascular Signal Processing Group
% Dresden, Germany, 2014
%
% daniel.wedekind@mailbox.tu-dresden.de
%
% Vers. 1.0:
% Last updated : 24-01-2014
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


filt_dat = zeros(size(data,1),size(data,2));
mov_dat = zeros(size(data,1),size(data,2));
diff_dat = zeros(size(data,1),size(data,2));

for index = 1:size(data,1)

    Nfir = 100;
    %bandstop around 50 Hz
    powerline = 50; %Hz
    BS = fir1(Nfir,[(powerline-1) (powerline+1)].*2/fs,'stop');
    BSdata = filtfilt(BS,1,data(index,:));
    %cascaded bandpass 5-11 Hz
    LP = fir1(Nfir,[11].*2/fs,'low');
    LPdata = filtfilt(LP,1,BSdata);
    HP = fir1(Nfir,[5].*2/fs,'high');
    HPdata = filtfilt(HP,1,LPdata);
    %smooth differentatition
    SD = [-(1/6),-(1/6),0,(1/6),(1/6)];
    SDdata = filtfilt(SD,1,HPdata);
    %squaring
    SQdata = SDdata.^2;
    %moving window integration of 150ms
    windowL = QRSmaxwidth+1;
    MOV = ones (1,windowL)/windowL;
    MOVINT = filtfilt(MOV,1,SQdata);
    
    filt_dat(index,:) = HPdata./std(HPdata);
    diff_dat(index,:) = SDdata./std(SDdata);
    mov_dat(index,:) = MOVINT./std(MOVINT);
    
end

return
