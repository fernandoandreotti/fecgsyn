function [ output_param ] = add_default_params( param , varargin)
%add_default_params Adds default parameter values
%   
%
% Mohsan Alvi (mohsan.alvi@eng.ox.ac.uk) - July 2014

% deal with optional inputs
if nargin < 2;    THR = 0.2;    else    THR = varargin{1};      end
if nargin < 3;    mVCG = 5;     else    mVCG = varargin{2};     end
if nargin < 4;    fVCG = 4;     else    fVCG = varargin{3};     end
if nargin < 5;    debug = 11;   else    debug = varargin{4};    end
if nargin < 6;    CH_CANC = 5;  else    CH_CANC = varargin{5};  end
if nargin < 7;    POS_DEV = 0;  else    POS_DEV = varargin{6};  end


% add default values if missing
% == default parameters
if ~any(strcmp('mheart',fieldnames(param))); param.mheart = [2*pi/3 0.2 0.4]; end;
if ~any(strcmp('fheart',fieldnames(param))); param.fheart{1} = [-pi/10 0.35 -0.3]; end;
if ~any(strcmp('elpos',fieldnames(param))); x = pi/12*[3 4 5 6 7 8 9 10]' -pi/2;     % 32 abdominal channels 
    y = .5*ones(8,1); xy = repmat([x y],4,1); z = repmat([-.1 -.2 -.3 -.4],8,1); z = reshape(z,32,1);
    abdmleads = [xy z]; refs = [-pi/4 0.5 0.4;(5/6-.5)*pi 0.5 0.4];  % + 2 reference leads
    param.elpos = [abdmleads;refs]; end   
if ~any(strcmp('refpos',fieldnames(param))); param.refpos = [pi 0.5 -0.3];end;
NB_FOETUSES = size(param.fheart,2); % number of foetuses figured out from the number of foetal heart locations entered
if ~any(strcmp('n',fieldnames(param))); param.n = 60000; end;
if ~any(strcmp('fs',fieldnames(param))); param.fs = 1000; end;
if ~any(strcmp('ntype',fieldnames(param))); param.ntype = ''; end;
if ~any(strcmp('noise_fct',fieldnames(param))); param.noise_fct(1:length(param.ntype)) = {1}; end;
if ~any(strcmp('SNRfm',fieldnames(param))); param.SNRfm = -9; end;
if ~any(strcmp('SNRmn',fieldnames(param))); param.SNRmn = 10; end;
if ~any(strcmp('mhr',fieldnames(param))); param.mhr = 90; end;
if ~any(strcmp('fhr',fieldnames(param))); param.fhr = repmat(150,NB_FOETUSES,1); end;
if ~any(strcmp('macc',fieldnames(param))); param.macc = 0; end;
if ~any(strcmp('facc',fieldnames(param))); param.facc = zeros(1,NB_FOETUSES); end;
if ~any(strcmp('mtypeacc',fieldnames(param))); param.mtypeacc = 'none'; end;
if ~any(strcmp('maccmean',fieldnames(param))); param.maccmean = 0; end;
if ~any(strcmp('maccstd',fieldnames(param))); param.maccstd = 1; end;
if ~any(strcmp('ftypeacc',fieldnames(param))); param.ftypeacc = arrayfun(@(x){sprintf('none',x)},1:NB_FOETUSES); end;
if ~any(strcmp('faccmean',fieldnames(param))); param.faccmean = repmat({0},1,NB_FOETUSES); end;
if ~any(strcmp('faccstd',fieldnames(param))); param.faccstd = repmat({1},1,NB_FOETUSES); end;
if ~any(strcmp('ftraj',fieldnames(param))); for cc=1:length(param.fhr); param.ftraj{cc} = 'none'; end; end;
if ~any(strcmp('fname',fieldnames(param))); param.fname = 'aecg'; end;
if ~any(strcmp('mres',fieldnames(param))); param.mres = 0; end;
if ~any(strcmp('fres',fieldnames(param))); param.fres = zeros(1,NB_FOETUSES); end;
if ~any(strcmp('mvcg',fieldnames(param))); param.mvcg = randi([1,9]); end;
if ~any(strcmp('fvcg',fieldnames(param))); param.fvcg = randi([1,9],NB_FOETUSES,1); end;
if ~any(strcmp('evcg',fieldnames(param))); param.evcg = randi([1,4]); end;
if ~any(strcmp('posdev',fieldnames(param))); param.posdev = 1; end;
if ~any(strcmp('mectb',fieldnames(param))); param.mectb = 0; end;
if ~any(strcmp('fectb',fieldnames(param))); param.fectb = 0; end;


output_param = param;

end

