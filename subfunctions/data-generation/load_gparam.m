function [gp,selvcg] = load_gparam(vcgmodel,type)
% load Gaussian paramters for vcg modelling. The parameters are derived 
% from the Physionet PTBDB and from the healthy control group on the first
% 20sec of the corresponding records. 
% PTBDB: http://www.physionet.org/physiobank/database/ptbdb/
%
% input:
%   vcgmodel: ID of Gaussian parameters to load
%   type:     type of ECG ('normal','ectopic')
%
% outputs:
%   gp: Gaussian parameters for ECG [cell]
%       gp{1} = tetai
%       gp{2} = alphai
%       gp{3} = bi 
%
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

vcgList = {'old1', 'old2', 's0273', 's0291', 's0302', 's0303', 's0306', 's0491', 's0533'};
selvcg = vcgList(vcgmodel);

if strcmp(type,'normal');
    % == normal ECG cycles case
    switch vcgmodel
            case 1
                tetai.x  = [-1.09  -0.83   -0.19     -.07  0 .06        0.22    1.2 1.42 1.68];
                alphai.x = [0.03   .08    -0.13    .85 1.11 .75     0.06   0.1  0.17 0.39];
                bi.x     = [0.0906    0.1057    0.0453    0.0378    0.0332    0.0302    0.0378    0.6040 0.3020    0.1812];

                tetai.y  = [-1.1  -0.9 -0.76       -0.11   -.01       0.065  0.8      1.58];
                alphai.y = [0.035 0.015 -0.019     0.32    .51     -0.32    0.04   0.08];
                bi.y     = [0.07  .07  0.04        0.055    0.037    0.0604  0.450  0.3];

                tetai.z  = [-1.1  -0.93 -0.7      -.4     -0.15    .095    1.05 1.25 1.55];
                alphai.z = [-0.03 -0.14 -0.035    .045     -0.4    .46    -.12 -.2 -.35];
                bi.z     = [.03  .12  .04         .4    .045       .05    .8 .4 .2];

            case 2    
                tetai.x  = [-0.7    -0.17    0       0.18     1.4];
                alphai.x = .1*[0.07     -0.11   1.3     0.07   0.275];
                bi.x     = [.1       .03     .045     0.02    0.3];

                tetai.y  = [-0.9     -0.08   0       0.05        1.3];
                alphai.y = [0.04     0.3     .45     -0.35       0.05];
                bi.y     = [.1       .05      .03    .04         .3];

                tetai.z  = [-0.8      -.3     -0.1        .06     1.35];
                alphai.z = .1*[-0.14    .03     -0.4        .46     -0.1];
                bi.z     = [.1       .4      .03         .03     .3];   

            case 3
                load('vcg_sets/GaussianParam_s0273lre_channel_x'); ParamX = Param;
                tetai.x  = ParamX(1,:);
                alphai.x = ParamX(2,:);
                bi.x     = ParamX(3,:);
                load('vcg_sets/GaussianParam_s0273lre_channel_y'); ParamY = Param;
                tetai.y  = ParamY(1,:);
                alphai.y = ParamY(2,:);
                bi.y     = ParamY(3,:);
                load('vcg_sets/GaussianParam_s0273lre_channel_z'); ParamZ = Param;
                tetai.z  = ParamZ(1,:);
                alphai.z = ParamZ(2,:);
                bi.z     = ParamZ(3,:);

            case 4
                load('vcg_sets/GaussianParam_s0291lre_channel_x'); ParamX = Param;
                tetai.x  = ParamX(1,:);
                alphai.x = ParamX(2,:);
                bi.x     = ParamX(3,:);
                load('vcg_sets/GaussianParam_s0291lre_channel_y'); ParamY = Param;
                tetai.y  = ParamY(1,:);
                alphai.y = ParamY(2,:);
                bi.y     = ParamY(3,:);
                load('vcg_sets/GaussianParam_s0291lre_channel_z'); ParamZ = Param;
                tetai.z  = ParamZ(1,:);
                alphai.z = ParamZ(2,:);
                bi.z     = ParamZ(3,:);           

            case 5
                load('vcg_sets/GaussianParam_s0302lre_channel_x'); ParamX = Param;
                tetai.x  = ParamX(1,:);
                alphai.x = ParamX(2,:);
                bi.x     = ParamX(3,:);
                load('vcg_sets/GaussianParam_s0302lre_channel_y'); ParamY = Param;
                tetai.y  = ParamY(1,:);
                alphai.y = ParamY(2,:);
                bi.y     = ParamY(3,:);
                load('vcg_sets/GaussianParam_s0302lre_channel_z'); ParamZ = Param;
                tetai.z  = ParamZ(1,:);
                alphai.z = ParamZ(2,:);
                bi.z     = ParamZ(3,:);     

            case 6
                load('vcg_sets/GaussianParam_s0303lre_channel_x'); ParamX = Param;
                tetai.x  = ParamX(1,:);
                alphai.x = ParamX(2,:);
                bi.x     = ParamX(3,:);
                load('vcg_sets/GaussianParam_s0303lre_channel_y'); ParamY = Param;
                tetai.y  = ParamY(1,:);
                alphai.y = ParamY(2,:);
                bi.y     = ParamY(3,:);
                load('vcg_sets/GaussianParam_s0303lre_channel_z'); ParamZ = Param;
                tetai.z  = ParamZ(1,:);
                alphai.z = ParamZ(2,:);
                bi.z     = ParamZ(3,:);    

            case 7
                load('vcg_sets/GaussianParam_s0306lre_channel_x'); ParamX = Param;
                tetai.x  = ParamX(1,:);
                alphai.x = ParamX(2,:);
                bi.x     = ParamX(3,:);
                load('vcg_sets/GaussianParam_s0306lre_channel_y'); ParamY = Param;
                tetai.y  = ParamY(1,:);
                alphai.y = ParamY(2,:);
                bi.y     = ParamY(3,:);
                load('vcg_sets/GaussianParam_s0306lre_channel_z'); ParamZ = Param;
                tetai.z  = ParamZ(1,:);
                alphai.z = ParamZ(2,:);
                bi.z     = ParamZ(3,:);             

            case 8
                load('vcg_sets/GaussianParam_s0491_re_channel_x'); ParamX = Param;
                tetai.x  = ParamX(1,:);
                alphai.x = ParamX(2,:);
                bi.x     = ParamX(3,:);
                load('vcg_sets/GaussianParam_s0491_re_channel_y'); ParamY = Param;
                tetai.y  = ParamY(1,:);
                alphai.y = ParamY(2,:);
                bi.y     = ParamY(3,:);
                load('vcg_sets/GaussianParam_s0491_re_channel_z'); ParamZ = Param;
                tetai.z  = ParamZ(1,:);
                alphai.z = ParamZ(2,:);
                bi.z     = ParamZ(3,:);

            case 9
                load('vcg_sets/GaussianParam_s0533_re_channel_x'); ParamX = Param;
                tetai.x  = ParamX(1,:);
                alphai.x = ParamX(2,:);
                bi.x     = ParamX(3,:);
                load('vcg_sets/GaussianParam_s0533_re_channel_y'); ParamY = Param;
                tetai.y  = ParamY(1,:);
                alphai.y = ParamY(2,:);
                bi.y     = ParamY(3,:);
                load('vcg_sets/GaussianParam_s0533_re_channel_z'); ParamZ = Param;
                tetai.z  = ParamZ(1,:);
                alphai.z = ParamZ(2,:);
                bi.z     = ParamZ(3,:);          
    end


elseif strcmp(type,'ectopic')
    % == Ectopic ad-hoc
    %switch vcgmodel
        %case 1
            tetai.x  = [-1.09,  -0.83,   -0.19,     -.07,  0, .06,        0.22,    1.2, 1.42, 1.68, 2.9];
            alphai.x = [0.03,   .08,    -0.13,    .65,              .70, .01,     0.06,   -0.05,  -0.17, -0.39, .03];
            bi.x     = [0.0906,    0.1057,    0.0453,    0.3378,    .5,    .2,    .5,    0.1040,  0.1020,   0.1812, .5];

            tetai.y  = [-1.1,  -0.9,   -0.76,     -.07,  0, .06,        0.22,    1.2, 1.42, 1.68, 2.9];
            alphai.y = [0.035,   .015,    -0.13,    .35,      .55 , .06,     0.06,   -0.05,  -0.17, -0.39, .014];
            bi.y     = [0.0906,    0.1057,    0.0453,    0.3378,    .5,    .2,    .5,    0.1040,  0.1020,   0.1812, .5];

            tetai.z  = [-1.09,  -0.83,   -0.19,     -.07,  0, .06,        0.22,    1.2, 1.42, 1.68, 2.9];
            alphai.z = [-0.03,   -.15,    -0.013,    .065,              -.70, .01,     -0.06,   0.05,  0.17, 0.39, .03];
            bi.z     = [0.0906,    0.1057,    0.0453,    0.3378,    .5,    .2,    .5,    0.1040,  0.1020,   0.1812, .5];
%         case 2
%             load('vcg_sets/EctopicBeatGaussians');
% 
%             alphai.x = OptimumParametersI2(1:7);
%             bi.x = OptimumParametersI2(8:14);
%             tetai.x = OptimumParametersI2(15:21);
% 
%             alphai.y = OptimumParametersI2(22:28);
%             bi.y = OptimumParametersI2(29:35);
%             tetai.y  = OptimumParametersI2(36:42);
% 
%             alphai.z = OptimumParametersI2(43:49);
%             bi.z = OptimumParametersI2(50:56);
%             tetai.z = OptimumParametersI2(57:63);    
%         case 3
%             load('vcg_sets/EctopicBeatGaussians');
% 
%             alphai.x = OptimumParametersI3(1:7);
%             bi.x = OptimumParametersI3(8:14);
%             tetai.x = OptimumParametersI3(15:21);
% 
%             alphai.y = OptimumParametersI3(22:28);
%             bi.y= OptimumParametersI3(29:35);
%             tetai.y  = OptimumParametersI3(36:42);
% 
%             alphai.z = OptimumParametersI3(43:49);
%             bi.z = OptimumParametersI3(50:56);
%             tetai.z = OptimumParametersI3(57:63);          
%         case 4
%             load('vcg_sets/EctopicBeatGaussians');
% 
%             alphai.x = OptimumParametersI35(1:7);
%             bi.x  = OptimumParametersI35(8:14);
%             tetai.x = OptimumParametersI35(15:21);
% 
%             alphai.y = OptimumParametersI35(22:28);
%             bi.y = OptimumParametersI35(29:35);
%             tetai.y = OptimumParametersI35(36:42);
% 
%             alphai.z = OptimumParametersI35(43:49);
%             bi.z = OptimumParametersI35(50:56);
%             tetai.z = OptimumParametersI35(57:63);  
%     end
else
   error('Add more cases if you want!') 
end

% == encapsule in structure
gp = {tetai, alphai, bi};

end






