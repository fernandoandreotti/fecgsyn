% FECGSYN_UI.
% FECGSYN Graphical user interface.
%
% FECGSYN is fruit of the collaboration between the Department of Engineering 
% Science, University of Oxford (DES-OX) and the Institute of Biomedical Engineering, 
% TU Dresden (IBMT-TUD). The authors are Joachim Behar (DES-OX), Fernando Andreotti 
% (IBMT-TUD), Julien Oster (DES-OX), Sebastian Zaunseder (IBMT-TUD) and 
% Gari Clifford (DES-OX). 
%
% The present user interface was contributed by Mohsan Alvi (DES-OX) under
% the supervision of Joachim Behar (DES-OX) and Fernando Andreotti (IBMT-TUD).
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

function varargout = FECGSYN_UI(varargin)
%% Initialization tasks

% Add files in parent folder and its subfolders to path
addpath('GUI');
addpath('noise_sources');
addpath('subfunctions');
addpath('vcg_sets');

mInputArgs = varargin;  % Command line arguments when invoking
                        % the GUI
mOutputArgs = {};       % Variable for storing output when GUI
                        % returns

% GUI mode: 1 - normal mode, 0 - debug mode using pregenerated 'out' file
gui_mode = 1;

% Create default font-size
if ismac
    fontSize = 14;
    fontSizeHelp = 11;
    fontSizeTitle = 24;
else
    fontSize = 13;
    fontSizeHelp = 8;
    fontSizeTitle = 16;
end
% Create constants for the calculations
THR = 0.2; % threshold of QRS detector
mVCG = 5; % choose mother VCG (if empty then the simulator randomly choose one within the set of available VCGs)
fVCG = 4; % choose foetus VCG (ibid)
debug = 11; % debug level. 11 corresponds to the gui
CH_CANC = 10; % max channel onto which to perform MECG cancellation
POS_DEV = 0; % slight deviation from default hearts and electrodes positions 
             % (0: hard coded values, 1: random deviation and phase initialisation)

% Initialise the "param" structures for each scenario
param_struct = FECGSYN_UI_create_param_structs(THR, mVCG, fVCG, debug, CH_CANC, POS_DEV);

% make 'out' global for import/export
out = struct;

% busy flag
busy_flag = 0;

% global variables holding the selected foetus and noise source so that the
% correct values can be saved when the user changes the selected object in
% the custom view list
selected_foetus = 1;
selected_ns = -1;

% Create master data for the default scenarios
num_scen = length(param_struct);
dscen_labels = cell(num_scen,1);
dscen_IDs = cell(num_scen,1);
for j = 1:num_scen
    dscen_labels{j} = param_struct{j}.title;
    dscen_IDs{j} = j;
end

% index of the popup menu
dscen_choice = 1; 

% The currently displayed plot's handle
active_plot_handle = [];

% Contains the list of plots that can be displayed
plot_titles = {};
plot_titles_children = {};

% All figure handles
f_handles = [];
electrode_plot_handle = []; % Used in the custom view's electrode preview 
displayed_electrode_plot_handle = [];

% Holds the "page" of elpos values being displayed. Each "page" shows up to
% 6 rows of the elpos matrix (out of a total of 34)
geo_num_labels = 3:34;  % choices in the popup menu choosing number of electrodes
elpos_page = 1;
elpos_page_size = 6; % number of electrodes shown on one page
elpos_displayed_idx = 1:elpos_page_size;

% Allowed choices for the dropdown variables in the custom view
ftypeacc_strings =   {'none' ...
                    , 'mexhat' ...
                    , 'gauss' ...
                    , 'flattop' ...
                    };

ftraj_strings =  {'none' ...
                , 'linear' ...
                , 'helix' ...
                , 'spline' ...
                , 'spiral' };

mtypeacc_strings =   {'none' ...
                    , 'mexhat' ...
                    , 'gauss' ...
                    , 'flattop' };
                    
ntype_strings =  {'MA' ...
                , 'EM' ...
                , 'BW' }; % list found in add_noisedipole.m
            
% noise_fct_strings =  {'1' ...
%                     , '2' ...
%                     , '3' ...
%                     , '4' ...
%                     , '5' };
                    
                    
                    
% Title stuff
gui_title = 'fecgsyn GUI';

% initialising the long texts
gui_text = {'Fecgsyn is a graphical user interface for the simulation of maternal-foetal activity mixtures on abdominal electrocardiogram recordings.'};

general_help_text = sprintf('Select a preset scenario from the dropdown menu, or create a custom simulation by pressing the custom button. Once a selection is made, click the run button to compute. The export button allows you to save a .mat file for results, and the import button opens a previously saved simulation.');

about_text = sprintf('FECGSYN is fruit of the collaboration between the Department of Engineering Science, University of Oxford (DES-OX) and the Institute of Biomedical Engineering, TU Dresden (IBMT-TUD). The authors are Joachim Behar (DES-OX), Fernando Andreotti (IBMT-TUD), Julien Oster (DES-OX), Sebastian Zaunseder (IBMT-TUD) and Gari Clifford (DES-OX). \n\nThe present user interface was contributed by Mohsan Alvi (DES-OX) under the supervision of Joachim Behar (DES-OX) and Fernando Andreotti (IBMT-TUD). \n\nThis program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details. \nYou should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>. \n \n Referencing this work \n\nfecgsyn simulator: \nBehar Joachim, Andreotti Fernando, Zaunseder Sebastian, Li Qiao, Oster Julien, Clifford Gari D. An ECG simulator for generating maternal-foetal activity mixtures on abdominal ECG recordings. Physiological Measurement.35 1537-1550. 2014. \n\nfecgsyn user interface software: \nAlvi Mohsan, Andreotti Fernando, Oster Julien, Sebastian Zaunseder , Clifford Gari D., Behar Joachim. A Graphical User Interface for the Simulation of Maternal-Foetal Activity Mixtures on Abdominal Electrocardiogram Recordings. Computing in Cardiology 2014 (accepted)');

custom_general_help_text = sprintf('General help \n \nn - number of samples \n \nfs - sampling frequency');

custom_fetal_help_text = sprintf('foetal help \n \nfheart - foetus heart origin - actual location will be picked randomly around it (default [-pi/10 0.4 -0.3])\n \nfhr - foetal heart rate [bpm] \n \nfacc - foetal acceleration in heart rate\n \nftypeacc - foetal acceleration type \n \nfectb - add ectopic beats to foetus [bool] \n \nfres - respiratory frequency of foetus [Hz] \n \nftraj - trajectory given to foetus heart \n \nfvcg - foetus vcg chosen (1-9)');

custom_noise_help_text = sprintf('Noise Help \n \nSNR_fm - Signal to noise ratio FECG/MECG (default -10) \n \nSNR_mn - Signal to noise ratio (MECG+FECG)/Noise (default 6) \n \nNoise source: \n \nntype - Noise type (default MA) \n \nnoise_fct - function of modulating noise (each noise may be modulated by a function of t, e.g sin(t))');

custom_mother_help_text = sprintf('Mother help \n \nmheart: maternal heart origin - actual location will be picked randomly around it. default [2*pi/3 0.4 0.4]) \n \nmhr - mother heart rate [bpm] \n \nmacc - maternal acceleration in heart rate \n \nmtypeacc - maternal acceleration type \n \nmectb - add ectopic beats to mother [bool] \n \nmres - respiratory frequency of mother [Hz] \n \nevcg - ectopic beat parameters (1-4)');

custom_controls_help_text = sprintf('Run - run fecgsyn \nEdit - edit a preset scenario \nBack - back to front face');

custom_scen_help_text = sprintf('Custom scenarios help \n \nDefault scenarios are provided for a variety of cases. These aim to illustrate the effect of different types of non-stationarities. If you wish to customize these, please press the edit button in the panel above.');

custom_geo_help_text = sprintf('Geometry help \n \nEnter electrode positions in spherical coordinates: \n \nth - rotation angle about z-axis \n \nz - height along z-axis \n \nThe preview button opens a graphical representation of the electrode positions. Click refresh to preview any changes made to electrode positions.');

%% Create handle to the GUI's main container
fh = figure('Name', 'fecgsyn GUI' ...      % Set title
          , 'NumberTitle', 'off' ...    % Hide "Figure 1:" from title
          , 'Position', [200 100 900 650] ...    % Position & size (x,y)
          , 'Resize', 'off' ...         % Disable resizing of the GUI
          , 'Visible', 'off' ...        % Start off with the GUI hidden
          , 'Tag', 'FECG_GUI' ...
          );
                
movegui(fh,'northwest');
      
%% Construct the main window components
fh_main = uipanel('Parent', fh ...
                    , 'Position', [0 0 1.0 1.0] ... 
                    , 'Visible', 'on');

                
% Main plot axis
axes_plots = uipanel('Parent', fh_main ...
                    , 'Position', [.02 .25 .96 .74] ...
                    );                
                                    
panel_about = uipanel('Parent', fh_main ...
                    , 'Position', [.02 .25 .96 .67] ...
                    , 'Visible', 'off');
                
list_plots = uicontrol(fh_main ...
                       , 'style', 'list' ...
                       , 'unit', 'pix' ...
                       , 'position', [43 20 245 120] ...
                       , 'min', 0 ...
                       , 'max', 1  ...
                       , 'Callback', @cb_list_plots ...
                       , 'string', plot_titles ...
                       , 'tooltip','Choose a plot to be displayed' ...
                       );

list_plots_children = uicontrol(fh_main ...
                       , 'style', 'list' ...
                       , 'unit', 'pix' ...
                       , 'position', [296 20 245 120] ...
                       , 'min', 0 ...
                       , 'max', 1  ...
                       , 'Callback', @cb_list_plots_children ...
                       , 'string', plot_titles_children ...
                       , 'tooltip','Choose a plot to be displayed' ...
                       );
                   
label_gui_busy = uicontrol(axes_plots,'Style','text'...
                               ,'String','Busy ...'...
                               ,'Visible', 'off' ...
                               ,'fontsize',20 ...                  
                               ,'Position',[350, 150, 100, 100]); 

% Default Scenario and Run
panel_defaults = uipanel('Parent', fh_main ...
                           , 'Position', [.65 .03 .3 .2] ...
                           , 'Title', 'Default Scenarios' ...
                           );
panel_defaults_help = uipanel('Parent', fh_main ...
                               , 'Position', [.35 .025 .6 .2] ...
                               , 'Title', get(panel_defaults, 'Title') ...
                               , 'Visible', 'off');    
                           
    % Pop-up menu to choose the scenario
    popup_default_scenario = uicontrol(panel_defaults,'Style','popupmenu',...
                    'String',dscen_labels, ...
                    'Value',1,'Position',[20 85 145 20], ...
                    'Callback',@cb_default_scenario_popup_menu, ...
                    'tooltip','Choose which scenario should be computed.');                
                
    % Run button 
    bt_run = uicontrol(panel_defaults...
                    ,'Style','pushbutton'...
                    ,'String','Run'...
                    ,'Position',[180 80 70 30]...
                    ,'Callback',@cb_run_button ...
                    ,'tooltip','Start the computation of the chosen scenario.' ...
            );

    % Import button 
    bt_import = uicontrol(panel_defaults...
                    ,'Style','pushbutton'...
                    ,'String','Import'...
                    ,'Position',[20 45 70 30]...
                    ,'Callback',@cb_import_button ...
                    ,'tooltip','Start the computation of the chosen scenario.' ...
            );
        
    % Export button 
    bt_export = uicontrol(panel_defaults...
                    ,'Style','pushbutton'...
                    ,'String','Export'...
                    ,'Position',[95 45 70 30]...
                    ,'Callback',@cb_export_button ...
                    ,'tooltip','Start the computation of the chosen scenario.' ...
            );
        
    % Custom view button    
    bt_open_custom = uicontrol(panel_defaults...
                    ,'Style','pushbutton'...
                    ,'String','Customize'...
                    ,'Position',[180 45 70 30]...
                    ,'Callback',@cb_open_custom ...
                    ,'tooltip','Customize the scenario parameters' ...
            );        

    % Exit button 
    bt_exit = uicontrol(panel_defaults...
                    ,'Style','pushbutton'...
                    ,'String','Exit'...
                    ,'Position',[180 10 70 30]...
                    ,'Callback','close all' ...
                    ,'tooltip','Close the GUI' ...
            );

    % Help button
    bt_defaults_help = uicontrol(panel_defaults...
                    ,'Style','pushbutton'...
                    ,'String','?'...
                    ,'Fontweight', 'bold' ...
                    ,'Position',[251 104 15 15]...
                    ,'Callback',@cb_bt_run_help ...
                    ,'tooltip', 'Help');
    % Help/Back button
    bt_defaults_help_x = uicontrol(panel_defaults_help...
                    ,'Style','pushbutton'...
                    ,'String','X'...
                    ,'Fontweight', 'bold' ...
                    ,'Position',[520 104 15 15]...
                    ,'Callback',@cb_bt_run_help_x ...
                    ,'tooltip', 'Help');
    
    % Help text
    label_defaults_help = uicontrol(panel_defaults_help,'Style','text'...
                               ,'String',general_help_text...
                               ,'fontsize',fontSizeHelp ...        
                               ,'HorizontalAlignment','left'...  
                               ,'Visible','on' ...        
                               ,'Position',[0 10 480 100]);     
                
    % About button - opens a message box with the about text
    bt_about = uicontrol(panel_defaults...
                    ,'Style','pushbutton'...
                    ,'String','About'...
                    ,'Position',[95 10 70 30]...
                    ,'Callback',@cb_bt_about ...
                    ,'tooltip','About' ...
            );
        

        
        
%% About page  

label_about_title = uicontrol(panel_about,'Style','text'...
                           ,'String',gui_title...
                           ,'fontsize',fontSizeTitle ...        
                           ,'HorizontalAlignment','left'...  
                           ,'Visible','on' ...        
                           ,'Position',[35 380 200 30]);
        
    label_about = uicontrol(panel_about,'Style','text'...
                           ,'String',about_text...
                           ,'fontsize',fontSizeHelp ...        
                           ,'HorizontalAlignment','left'...  
                           ,'Visible','on' ...        
                           ,'Position',[320 20 470 390]);

    bt_about_close = uicontrol(panel_about...
                    ,'Style','pushbutton'...
                    ,'String','Back'...
                    ,'Position',[680 10 70 30]...
                    ,'Callback',@cb_bt_about_close ...
                    ,'tooltip','Close About' ...
            );
        
    logo1_axes = axes('Parent', panel_about ...
                    , 'Position', [0.03 0.46 0.2 0.5] ...
                    );
    axes(logo1_axes)
    imshow('FrenchFlag.png')

    logo2_axes = axes('Parent', panel_about ...
                    , 'Position', [0.03 0.18 0.3 0.55] ...
                    );
    axes(logo2_axes)
    imshow('oxlogo_big.png') 
    
    logo3_axes = axes('Parent', panel_about ...
                    , 'Position', [0.03 0 0.3 0.5] ...
                    );
    axes(logo3_axes)
    imshow('tu_logo.jpg')
                           
%% Construct the custom window components
fh_custom = uipanel('Parent', fh ...
                    , 'Position', [0 0 1.0 1.0] ...
                    , 'Visible', 'off');

% Panels for parameter groups (naming convention: panel_'panel name'_params)

    panel_general_params = uipanel('Parent', fh_custom ...
                               , 'Position', [0.02 .7 .3 .2] ...
                               , 'Title', 'General Params');

    panel_fetal_params = uipanel('Parent', fh_custom ...
                               , 'Position', [.02 .1 .3 .6] ...
                               , 'Title', 'Foetal Params');

    panel_noise_params = uipanel('Parent', fh_custom ...
                               , 'Position', [.34 .5 .3 .4] ...
                               , 'Title', 'Noise Params');

    panel_mother_params = uipanel('Parent', fh_custom ...
                               , 'Position', [0.34 .1 .3 .4] ...
                               , 'Title', 'Mother Params');
                           
    panel_elpos_preview = uipanel('Parent', fh_custom ...
                               , 'Position', [0.02 .1 .62 .8] ...
                               , 'Title', 'Electrode Positions'...
                               , 'Visible', 'off');
                           
    panel_custom_controls = uipanel('Parent', fh_custom ...
                               , 'Position', [.66 .80 .3 .1] ...
                               , 'Title', 'Controls');
                           
    panel_custom_scenarios = uipanel('Parent', fh_custom ...
                               , 'Position', [.66 .5 .3 .3] ...
                               , 'Title', 'Load Scenario Parameters');

    panel_geometry = uipanel('Parent', fh_custom ...
                               , 'Position', [.66 .1 .3 .4] ...
                               , 'Title', 'Geometry Params');                       

    panel_geometry_no_edit = uipanel('Parent', fh_custom ...
                               , 'Position', get(panel_geometry, 'Position') ...
                               , 'Title', 'Geometry Params' ...
                               , 'Visible','off');
                           
    % Custom view Help boxes
    panel_general_params_help = uipanel('Parent', fh_custom ...
                               , 'Position', get(panel_general_params, 'Position') ...
                               , 'Title', get(panel_general_params, 'Title') ...
                               , 'Visible', 'off');

    panel_fetal_params_help = uipanel('Parent', fh_custom ...
                               , 'Position', get(panel_fetal_params, 'Position') ...
                               , 'Title', get(panel_fetal_params, 'Title') ...
                               , 'Visible', 'off');

    panel_noise_params_help = uipanel('Parent', fh_custom ...
                               , 'Position', get(panel_noise_params, 'Position') ...
                               , 'Title', get(panel_noise_params, 'Title') ...
                               , 'Visible', 'off');

    panel_mother_params_help = uipanel('Parent', fh_custom ...
                               , 'Position', get(panel_mother_params, 'Position') ...
                               , 'Title', get(panel_mother_params, 'Title') ...
                               , 'Visible', 'off');
                           
    panel_custom_controls_help = uipanel('Parent', fh_custom ...
                               , 'Position', get(panel_custom_controls, 'Position') ...
                               , 'Title', get(panel_custom_controls, 'Title') ...
                               , 'Visible', 'off');
                           
    panel_custom_scenarios_help = uipanel('Parent', fh_custom ...
                               , 'Position', get(panel_custom_scenarios, 'Position') ...
                               , 'Title', get(panel_custom_scenarios, 'Title') ...
                               , 'Visible', 'off');

    panel_geometry_help = uipanel('Parent', fh_custom ...
                               , 'Position', get(panel_geometry, 'Position') ...
                               , 'Title', get(panel_geometry, 'Title') ...
                               , 'Visible', 'off');        
                           
                           
    % Main controls
    bt_run_custom = uicontrol(panel_custom_controls...
                    ,'Style','pushbutton'...
                    ,'String','Run'...
                    ,'Position',[15 7 70 40]...
                    ,'Callback',@cb_bt_run_custom ...
            );  
        
    bt_save_edit_custom = uicontrol(panel_custom_controls...
                    ,'Style','pushbutton'...
                    ,'String','Edit'...
                    ,'Position',[95 7 70 40]...
                    ,'Callback',@cb_save_edit_custom ...
            );

    bt_back_custom = uicontrol(panel_custom_controls...
                    ,'Style','pushbutton'...
                    ,'String','Back'...
                    ,'Position',[175 7 70 40]...
                    ,'Callback',@cb_back_custom ...
            );    
% Input fields (naming convention: input_'name of input')
% The positions of input boxes is set starting at y_input from the bottom,
% x_input from the left, at distance y_input_diff to each other vertically. Each box is 
% x_input_width by y_input_height in size.
y_input = 20; y_input_diff = 5; y_input_height = 20; x_input = 120; x_input_width = 130;
x_offset = 5;
% buffer size for spacing between nearby components
buf = 3;
% size of text fields if there are two or three in a row
x_input_width_2 = (x_input_width - buf)/2;
x_input_width_3 = (x_input_width - 2*buf)/3;

% Mother parameters
    input_mother_1_1 = uicontrol(panel_mother_params, 'Style', 'edit' ...
                               ,'Position',[x_input, y_input+(y_input_diff+y_input_height)*7, x_input_width_3, y_input_height] ...
                               , 'KeyPressFcn',@cb_enter_key);
    input_mother_1_2 = uicontrol(panel_mother_params, 'Style', 'edit' ...
                               ,'Position',[x_input+(x_input_width_3+buf), y_input+(y_input_diff+y_input_height)*7, x_input_width_3, y_input_height] ...
                               , 'KeyPressFcn',@cb_enter_key);
    input_mother_1_3 = uicontrol(panel_mother_params, 'Style', 'edit' ...
                               ,'Position',[x_input+2*(x_input_width_3+buf), y_input+(y_input_diff+y_input_height)*7, x_input_width_3, y_input_height] ...
                               , 'KeyPressFcn',@cb_enter_key);
    label_mother_1 = uicontrol(panel_mother_params,'Style','text'...
                               ,'String','mheart'...
                               ,'fontsize',fontSize ...
                               ,'HorizontalAlignment','left'...
                               ,'Position',[x_offset, y_input+(y_input_diff+y_input_height)*7-2, x_input_width-20, y_input_height]);

    input_mother_2 = uicontrol(panel_mother_params, 'Style', 'edit' ...
                               ,'Position',[x_input, y_input+(y_input_diff+y_input_height)*6, x_input_width, y_input_height] ...
                               , 'KeyPressFcn',@cb_enter_key);
    label_mother_2 = uicontrol(panel_mother_params,'Style','text'...
                               ,'String','mhr'...
                               ,'fontsize',fontSize ...
                               ,'HorizontalAlignment','left'...
                               ,'Position',[x_offset, y_input+(y_input_diff+y_input_height)*6-2, x_input_width-20, y_input_height]);

    input_mother_3 = uicontrol(panel_mother_params, 'Style', 'edit' ...
                               ,'Position',[x_input, y_input+(y_input_diff+y_input_height)*5, x_input_width, y_input_height] ...
                               , 'KeyPressFcn',@cb_enter_key);  
    label_mother_3 = uicontrol(panel_mother_params,'Style','text'...
                               ,'String','macc'...
                               ,'fontsize',fontSize ...
                               ,'HorizontalAlignment','left'...
                               ,'Position',[x_offset, y_input+(y_input_diff+y_input_height)*5-2, x_input_width-20, y_input_height]);  

    %input_mother_4 = uicontrol(panel_mother_params, 'Style', 'edit' ...
%                                ,'Position',[x_input, y_input+(y_input_diff+y_input_height)*4, x_input_width, y_input_height] ...
%                                , 'KeyPressFcn',@cb_enter_key); 
    popup_mother_mtypeacc = uicontrol(panel_mother_params,'Style','popupmenu',...
                    'String',mtypeacc_strings, ...
                    'Value',1, ...
                    'Position',[x_input, y_input+(y_input_diff+y_input_height)*4, x_input_width, y_input_height], ...
                    'Callback',@cb_popup_mother_mtypeacc);
    label_mother_4 = uicontrol(panel_mother_params,'Style','text'...
                               ,'String','mtypeacc'...
                               ,'fontsize',fontSize ...   
                               ,'HorizontalAlignment','left'...                        
                               ,'Position',[x_offset, y_input+(y_input_diff+y_input_height)*4-2, x_input_width-20, y_input_height]); 

%     input_mother_5 = uicontrol(panel_mother_params, 'Style', 'edit' ...
%                                ,'Position',[x_input, y_input+(y_input_diff+y_input_height)*3, x_input_width, y_input_height] ...
%                                , 'KeyPressFcn',@cb_enter_key);
    chk_mother_mecbt = uicontrol(panel_mother_params, 'Style', 'checkbox' ...
                               ,'String', '' ...
                               ,'Position',[x_input+4/3*x_input_width_3+buf, y_input+(y_input_diff+y_input_height)*3, x_input_width_3, y_input_height] ...
                               , 'Callback',@cb_chk_fetal_fecbt);
    label_mother_5 = uicontrol(panel_mother_params,'Style','text'...
                               ,'String','mectb'...
                               ,'fontsize',fontSize ...
                               ,'HorizontalAlignment','left'...
                               ,'Position',[x_offset, y_input+(y_input_diff+y_input_height)*3-2, x_input_width-20, y_input_height]);

    input_mother_6 = uicontrol(panel_mother_params, 'Style', 'edit' ...
                               ,'Position',[x_input, y_input+(y_input_diff+y_input_height)*2, x_input_width, y_input_height] ...
                               , 'KeyPressFcn',@cb_enter_key); 
    label_mother_6 = uicontrol(panel_mother_params,'Style','text'...
                               ,'String','mres'...
                               ,'fontsize',fontSize ...
                               ,'HorizontalAlignment','left'...
                               ,'Position',[x_offset, y_input+(y_input_diff+y_input_height)*2-2, x_input_width-20, y_input_height]);    

    input_mother_7 = uicontrol(panel_mother_params, 'Style', 'edit' ...
                               ,'Position',[x_input, y_input+(y_input_diff+y_input_height), x_input_width, y_input_height] ...
                               , 'KeyPressFcn',@cb_enter_key);
    label_mother_7 = uicontrol(panel_mother_params,'Style','text'...
                               ,'String','evcg'...
                               ,'fontsize',fontSize ...
                               ,'HorizontalAlignment','left'...
                               ,'Position',[x_offset, y_input+(y_input_diff+y_input_height)-2, x_input_width-20, y_input_height]);
                           

% General parameters
    input_general_1 = uicontrol(panel_general_params, 'Style', 'edit' ...
                               ,'Position',[x_input, y_input+(y_input_diff+y_input_height)*2, x_input_width, y_input_height] ...
                               , 'KeyPressFcn',@cb_enter_key );
    label_general_1 = uicontrol(panel_general_params,'Style','text'...
                               ,'String','n'...
                               ,'fontsize',fontSize ...      
                               ,'HorizontalAlignment','left'...            
                               ,'Position',[x_offset, y_input+(y_input_diff+y_input_height)*2-2, x_input_width-20, y_input_height]);
                           
    input_general_2 = uicontrol(panel_general_params, 'Style', 'edit' ...
                               ,'Position',[x_input, y_input+(y_input_diff+y_input_height)*1, x_input_width, y_input_height] ...
                               , 'KeyPressFcn',@cb_enter_key);
    label_general_2 = uicontrol(panel_general_params,'Style','text'...
                               ,'String','fs'...
                               ,'fontsize',fontSize ...     
                               ,'HorizontalAlignment','left'...             
                               ,'Position',[x_offset, y_input+(y_input_diff+y_input_height)*1-2, x_input_width-20, y_input_height]);
    
% Noise parameters
    input_noise_1 = uicontrol(panel_noise_params, 'Style', 'edit' ...
                               ,'Position',[x_input, y_input+(y_input_diff+y_input_height)*8, x_input_width, y_input_height] ...
                               , 'KeyPressFcn',@cb_enter_key);
    label_noise_1 = uicontrol(panel_noise_params,'Style','text'...
                               ,'String','SNR_fm'...
                               ,'fontsize',fontSize ...                
                               ,'HorizontalAlignment','left'...  
                               ,'Position',[x_offset, y_input+(y_input_diff+y_input_height)*8-2, x_input_width-20, y_input_height]); 
                           
    input_noise_2 = uicontrol(panel_noise_params, 'Style', 'edit' ...
                               ,'Position',[x_input, y_input+(y_input_diff+y_input_height)*7, x_input_width, y_input_height] ...
                               , 'KeyPressFcn',@cb_enter_key);
    label_noise_2 = uicontrol(panel_noise_params,'Style','text'...
                               ,'String','SNR_mn'...
                               ,'fontsize',fontSize ...                
                               ,'HorizontalAlignment','left'... 
                               ,'Position',[x_offset, y_input+(y_input_diff+y_input_height)*7-2, x_input_width-20, y_input_height]);                            
            
    label_noise = uicontrol(panel_noise_params,'Style','text'...
                               ,'String','Noise Source:'...
                               ,'fontsize',fontSize ...                
                               ,'Position',[x_offset, y_input+(y_input_diff+y_input_height)*6-2, x_input_width-20, y_input_height]);                             
                           
    %input_noise_3 = uicontrol(panel_noise_params, 'Style', 'edit' ...
%                                ,'Position',[x_input, y_input+(y_input_diff+y_input_height)*5, x_input_width, y_input_height] ...
%                                , 'KeyPressFcn',@cb_enter_key);
    popup_noise_ntype = uicontrol(panel_noise_params,'Style','popupmenu',...
                    'String',ntype_strings, ...
                    'Value',1, ...
                    'Position',[x_input, y_input+(y_input_diff+y_input_height)*5, x_input_width, y_input_height], ...
                    'Callback',@cb_popup_noise_ntype);
    
    label_noise_3 = uicontrol(panel_noise_params,'Style','text'...
                               ,'String','ntype'...
                               ,'fontsize',fontSize ...         
                               ,'HorizontalAlignment','left'...         
                               ,'Position',[x_offset, y_input+(y_input_diff+y_input_height)*5-2, x_input_width-20, y_input_height]);                            
                           

    input_noise_4 = uicontrol(panel_noise_params, 'Style', 'edit' ...
                                ,'Position',[x_input, y_input+(y_input_diff+y_input_height)*4, x_input_width, y_input_height] ...
                                , 'KeyPressFcn',@cb_enter_key); 
%     popup_noise_fct = uicontrol(panel_noise_params,'Style','popupmenu',...
%                     'String',noise_fct_strings, ...
%                     'Value',1, ...
%                     'Position',[x_input, y_input+(y_input_diff+y_input_height)*4, x_input_width, y_input_height], ...
%                     'Callback',@cb_popup_noise_fct);
    label_noise_4 = uicontrol(panel_noise_params,'Style','text'...
                               ,'String','noise_fct'...
                               ,'fontsize',fontSize ...      
                               ,'HorizontalAlignment','left'...            
                               ,'Position',[x_offset, y_input+(y_input_diff+y_input_height)*4-2, x_input_width-20, y_input_height]);                            
                           
                           
% Geometry params
    label_geo_num = uicontrol(panel_geometry,'Style','text'...
                               ,'String', 'Number of electrodes' ... 
                               ,'fontsize',fontSize ...
                               ,'Position',[0, y_input+(y_input_diff+y_input_height)*8, 180, y_input_height]);
    
    % Pop-up menu to choose the number of electrodes
    popup_geo_num = uicontrol(panel_geometry,'Style','popupmenu',...
                    'String',geo_num_labels ...
                    , 'Value',length(geo_num_labels) ...
                    , 'Position',[x_input+(x_input_width_2), y_input+(y_input_diff+y_input_height)*8, x_input_width_2, y_input_height] ...
                    , 'Callback',@cb_popup_geo_num);
    
    label_geo_th = uicontrol(panel_geometry,'Style','text'...
                               ,'String', 'q' ... 
                               , 'FontName','symbol' ...
                               ,'fontsize',fontSize ...
                               ,'Position',[x_input, y_input+(y_input_diff+y_input_height)*7, x_input_width_2, y_input_height]);
    label_geo_z = uicontrol(panel_geometry,'Style','text'...
                                   ,'String','z'...
                                   ,'fontsize',fontSize ...
                                   ,'Position',[x_input+(x_input_width_2+buf), y_input+(y_input_diff+y_input_height)*7, x_input_width_2, y_input_height]);
%     label_geo_z = uicontrol(panel_geometry,'Style','text'...
%                                    ,'String','z'...
%                                    ,'fontsize',fontSize ...
%                                    ,'Position',[x_input+2*(x_input_width_2+buf), y_input+(y_input_diff+y_input_height)*8, x_input_width_2, y_input_height]);

    input_geo_1_1 = uicontrol(panel_geometry, 'Style', 'edit' ...
                               ,'Position',[x_input, y_input+(y_input_diff+y_input_height)*6, x_input_width_2, y_input_height] ...
                               , 'KeyPressFcn',@cb_enter_key);
%     input_geo_1_2 = uicontrol(panel_geometry, 'Style', 'edit' ...
%                                ,'Position',[x_input+(x_input_width_3+buf), y_input+(y_input_diff+y_input_height)*7, x_input_width_3, y_input_height] ...
%                                , 'KeyPressFcn',@cb_enter_key);
    input_geo_1_3 = uicontrol(panel_geometry, 'Style', 'edit' ...
                               ,'Position',[x_input+(x_input_width_2+buf), y_input+(y_input_diff+y_input_height)*6, x_input_width_2, y_input_height] ...
                               , 'KeyPressFcn',@cb_enter_key);
    label_geo_1 = uicontrol(panel_geometry,'Style','text'...
                               ,'String','Electrode 1'...
                               ,'fontsize',fontSize ...
                               ,'HorizontalAlignment','left'...
                               ,'Position',[x_offset, y_input+(y_input_diff+y_input_height)*6+2, x_input_width-20, y_input_height]);

    input_geo_2_1 = uicontrol(panel_geometry, 'Style', 'edit' ...
                               ,'Position',[x_input, y_input+(y_input_diff+y_input_height)*5, x_input_width_2, y_input_height] ...
                               , 'KeyPressFcn',@cb_enter_key);
%     input_geo_2_2 = uicontrol(panel_geometry, 'Style', 'edit' ...
%                                ,'Position',[x_input+(x_input_width_3+buf), y_input+(y_input_diff+y_input_height)*6, x_input_width_3, y_input_height] ...
%                                , 'KeyPressFcn',@cb_enter_key);
    input_geo_2_3 = uicontrol(panel_geometry, 'Style', 'edit' ...
                               ,'Position',[x_input+(x_input_width_2+buf), y_input+(y_input_diff+y_input_height)*5, x_input_width_2, y_input_height] ...
                               , 'KeyPressFcn',@cb_enter_key);
    label_geo_2 = uicontrol(panel_geometry,'Style','text'...
                               ,'String','Electrode 2'...
                               ,'fontsize',fontSize ...
                               ,'HorizontalAlignment','left'...
                               ,'Position',[x_offset, y_input+(y_input_diff+y_input_height)*5+2, x_input_width-20, y_input_height]);
                           
    input_geo_3_1 = uicontrol(panel_geometry, 'Style', 'edit' ...
                               ,'Position',[x_input, y_input+(y_input_diff+y_input_height)*4, x_input_width_2, y_input_height] ...
                               , 'KeyPressFcn',@cb_enter_key);
%     input_geo_3_2 = uicontrol(panel_geometry, 'Style', 'edit' ...
%                                ,'Position',[x_input+(x_input_width_3+buf), y_input+(y_input_diff+y_input_height)*5, x_input_width_3, y_input_height] ...
%                                , 'KeyPressFcn',@cb_enter_key);
    input_geo_3_3 = uicontrol(panel_geometry, 'Style', 'edit' ...
                               ,'Position',[x_input+(x_input_width_2+buf), y_input+(y_input_diff+y_input_height)*4, x_input_width_2, y_input_height] ...
                               , 'KeyPressFcn',@cb_enter_key);
    label_geo_3 = uicontrol(panel_geometry,'Style','text'...
                               ,'String','Electrode 3'...
                               ,'fontsize',fontSize ...
                               ,'HorizontalAlignment','left'...
                               ,'Position',[x_offset, y_input+(y_input_diff+y_input_height)*4+2, x_input_width-20, y_input_height]);

    input_geo_4_1 = uicontrol(panel_geometry, 'Style', 'edit' ...
                               ,'Position',[x_input, y_input+(y_input_diff+y_input_height)*3, x_input_width_2, y_input_height] ...
                               , 'KeyPressFcn',@cb_enter_key);
%     input_geo_4_2 = uicontrol(panel_geometry, 'Style', 'edit' ...
%                                ,'Position',[x_input+(x_input_width_3+buf), y_input+(y_input_diff+y_input_height)*4, x_input_width_3, y_input_height] ...
%                                , 'KeyPressFcn',@cb_enter_key);
    input_geo_4_3 = uicontrol(panel_geometry, 'Style', 'edit' ...
                               ,'Position',[x_input+(x_input_width_2+buf), y_input+(y_input_diff+y_input_height)*3, x_input_width_2, y_input_height] ...
                               , 'KeyPressFcn',@cb_enter_key);
    label_geo_4 = uicontrol(panel_geometry,'Style','text'...
                               ,'String','Electrode 4'...
                               ,'fontsize',fontSize ...
                               ,'HorizontalAlignment','left'...
                               ,'Position',[x_offset, y_input+(y_input_diff+y_input_height)*3+2, x_input_width-20, y_input_height]);
                           
    input_geo_5_1 = uicontrol(panel_geometry, 'Style', 'edit' ...
                               ,'Position',[x_input, y_input+(y_input_diff+y_input_height)*2, x_input_width_2, y_input_height] ...
                               , 'KeyPressFcn',@cb_enter_key);
%     input_geo_5_2 = uicontrol(panel_geometry, 'Style', 'edit' ...
%                                ,'Position',[x_input+(x_input_width_3+buf), y_input+(y_input_diff+y_input_height)*3, x_input_width_3, y_input_height] ...
%                                , 'KeyPressFcn',@cb_enter_key);
    input_geo_5_3 = uicontrol(panel_geometry, 'Style', 'edit' ...
                               ,'Position',[x_input+(x_input_width_2+buf), y_input+(y_input_diff+y_input_height)*2, x_input_width_2, y_input_height] ...
                               , 'KeyPressFcn',@cb_enter_key);
    label_geo_5 = uicontrol(panel_geometry,'Style','text'...
                               ,'String','Electrode 5'...
                               ,'fontsize',fontSize ...
                               ,'HorizontalAlignment','left'...
                               ,'Position',[x_offset, y_input+(y_input_diff+y_input_height)*2+2, x_input_width-20, y_input_height]);

    input_geo_6_1 = uicontrol(panel_geometry, 'Style', 'edit' ...
                               ,'Position',[x_input, y_input+(y_input_diff+y_input_height)*1, x_input_width_2, y_input_height] ...
                               , 'KeyPressFcn',@cb_enter_key);
%     input_geo_6_2 = uicontrol(panel_geometry, 'Style', 'edit' ...
%                                ,'Position',[x_input+(x_input_width_3+buf), y_input+(y_input_diff+y_input_height)*2, x_input_width_3, y_input_height] ...
%                                , 'KeyPressFcn',@cb_enter_key);
    input_geo_6_3 = uicontrol(panel_geometry, 'Style', 'edit' ...
                               ,'Position',[x_input+(x_input_width_2+buf), y_input+(y_input_diff+y_input_height)*1, x_input_width_2, y_input_height] ...
                               , 'KeyPressFcn',@cb_enter_key);
    label_geo_6 = uicontrol(panel_geometry,'Style','text'...
                               ,'String','Electrode 6'...
                               ,'fontsize',fontSize ...
                               ,'HorizontalAlignment','left'...
                               ,'Position',[x_offset, y_input+(y_input_diff+y_input_height)*1+2, x_input_width-20, y_input_height]);
                           
%     input_geo_7_1 = uicontrol(panel_geometry, 'Style', 'edit' ...
%                                ,'Position',[x_input, y_input+(y_input_diff+y_input_height)*1, x_input_width_2, y_input_height] ...
%                                , 'KeyPressFcn',@cb_enter_key);
%     input_geo_7_2 = uicontrol(panel_geometry, 'Style', 'edit' ...
%                                ,'Position',[x_input+(x_input_width_3+buf), y_input+(y_input_diff+y_input_height)*1, x_input_width_3, y_input_height] ...
%                                , 'KeyPressFcn',@cb_enter_key);
%     input_geo_7_3 = uicontrol(panel_geometry, 'Style', 'edit' ...
%                                ,'Position',[x_input+(x_input_width_2+buf), y_input+(y_input_diff+y_input_height)*1, x_input_width_2, y_input_height] ...
%                                , 'KeyPressFcn',@cb_enter_key);
%     label_geo_7 = uicontrol(panel_geometry,'Style','text'...
%                                ,'String','Electrode 7'...
%                                ,'fontsize',fontSize ...
%                                ,'HorizontalAlignment','left'...
%                                ,'Position',[x_offset, y_input+(y_input_diff+y_input_height)*1+2, x_input_width-20, y_input_height]);
      
    bt_geo_previous = uicontrol(panel_geometry...
                    ,'Style','pushbutton'...
                    ,'String','<'...
                    ,'Position',[15 5 50 30]...
                    ,'Callback',@cb_bt_geo_previous ...
            );        
    bt_geo_preview = uicontrol(panel_geometry...
                    ,'Style','pushbutton'...
                    ,'String','Preview'...
                    ,'Position',[95 5 80 30]...
                    ,'Callback',@cb_bt_geo_preview ...
            ); 
        
    bt_geo_next = uicontrol(panel_geometry...
                    ,'Style','pushbutton'...
                    ,'String','>'...
                    ,'Position',[205 5 50 30]...
                    ,'Callback',@cb_bt_geo_next ...
            );  
    
    bt_geo_hide = uicontrol(panel_geometry...
                    ,'Style','pushbutton'...
                    ,'String','Hide'...
                    ,'Position',[75 5 55 30]...
                    ,'Visible', 'off' ...
                    ,'Callback',@cb_bt_geo_hide ...
            ); 
    
    bt_geo_refresh = uicontrol(panel_geometry...
                    ,'Style','pushbutton'...
                    ,'String','Refresh'...
                    ,'Position',[140 5 55 30]...
                    ,'Visible', 'off' ...
                    ,'Callback',@cb_bt_geo_refresh ...
            ); 
        
    geo_ui_elements = [ label_geo_1, input_geo_1_1, input_geo_1_3;
                        label_geo_2, input_geo_2_1, input_geo_2_3;
                        label_geo_3, input_geo_3_1, input_geo_3_3;
                        label_geo_4, input_geo_4_1, input_geo_4_3;
                        label_geo_5, input_geo_5_1, input_geo_5_3;
                        label_geo_6, input_geo_6_1, input_geo_6_3;
%                         label_geo_7, input_geo_7_1, input_geo_7_3;
                        ];
                    
    bt_geo_hide_no_edit = uicontrol(panel_geometry_no_edit...
                    ,'Style','pushbutton'...
                    ,'String','Hide'...
                    ,'Position',[95 5 80 30]...
                    ,'Visible', 'off' ...
                    ,'Callback',@cb_bt_geo_hide ...
            ); 
                    
    bt_geo_edit = uicontrol(panel_geometry_no_edit...
                    ,'Style','pushbutton'...
                    ,'String','Edit'...
                    ,'Position',[95 45 80 30]...
                    ,'Callback',@cb_save_edit_custom ...
            );
        
    bt_geo_preview_no_edit = uicontrol(panel_geometry_no_edit...
                    ,'Style','pushbutton'...
                    ,'String','Preview'...
                    ,'Position',[95 85 80 30]...
                    ,'Callback',@cb_bt_geo_preview ...
            ); 
                    
    % List of noise sources. These can be added and removed using buttons add foetus and remove foetus. Adding a noise source
    % will create an entry in this list. When the noise source is selected, its parameters should be visible.
    list_noise_sources = uicontrol(panel_noise_params,'style','list',...
                     'unit','pix',...
                     'position',[10 50 245 60],...
                     'min',0,'max',2,...
                     'fontsize',fontSize,...
                     'callback', @cb_list_noise_sources,...
                     'string',{'noise source 1', 'noise source 2'});

    % Add Noise source button
    bt_add_noise = uicontrol(panel_noise_params...
                    ,'Style','pushbutton'...
                    ,'String','Add Noise Source'...
                    ,'Position',[10 10 100 30]...
                    ,'Callback',@cb_add_noise ...
            );

    % Remove Noise source button
    bt_remove_noise = uicontrol(panel_noise_params...
                    ,'Style','pushbutton'...
                    ,'String','Remove Noise Source'...
                    ,'Position',[120 10 100 30]...
                    ,'Callback',@cb_remove_noise ...
        );   

    
% Fetal parameters
    input_fetal_1_1 = uicontrol(panel_fetal_params, 'Style', 'edit' ...
                               ,'Position',[x_input, y_input+(y_input_diff+y_input_height)*13, x_input_width_3, y_input_height] ...
                               , 'KeyPressFcn',@cb_enter_key);
    input_fetal_1_2 = uicontrol(panel_fetal_params, 'Style', 'edit' ...
                               ,'Position',[x_input+x_input_width_3+buf, y_input+(y_input_diff+y_input_height)*13, x_input_width_3, y_input_height] ...
                               , 'KeyPressFcn',@cb_enter_key);
    input_fetal_1_3 = uicontrol(panel_fetal_params, 'Style', 'edit' ...
                               ,'Position',[x_input+2*(x_input_width_3+buf), y_input+(y_input_diff+y_input_height)*13, x_input_width_3, y_input_height] ...
                               , 'KeyPressFcn',@cb_enter_key);
    label_fetal_1 = uicontrol(panel_fetal_params,'Style','text'...
                               ,'String','fheart'...
                               ,'fontsize',fontSize ...   
                               ,'HorizontalAlignment','left'...               
                               ,'Position',[x_offset, y_input+(y_input_diff+y_input_height)*13-2, x_input_width-20, y_input_height]); 
                           
    input_fetal_2 = uicontrol(panel_fetal_params, 'Style', 'edit' ...
                               ,'Position',[x_input, y_input+(y_input_diff+y_input_height)*12, x_input_width, y_input_height] ...
                               , 'KeyPressFcn',@cb_enter_key);
    label_fetal_2 = uicontrol(panel_fetal_params,'Style','text'...
                               ,'String','fhr'...
                               ,'fontsize',fontSize ...    
                               ,'HorizontalAlignment','left'...              
                               ,'Position',[x_offset, y_input+(y_input_diff+y_input_height)*12-2, x_input_width-20, y_input_height]); 
                           
    input_fetal_3 = uicontrol(panel_fetal_params, 'Style', 'edit' ...
                               ,'Position',[x_input, y_input+(y_input_diff+y_input_height)*11, x_input_width, y_input_height] ...
                               , 'KeyPressFcn',@cb_enter_key);  
    label_fetal_3 = uicontrol(panel_fetal_params,'Style','text'...
                               ,'String','facc'...
                               ,'fontsize',fontSize ...        
                               ,'HorizontalAlignment','left'...          
                               ,'Position',[x_offset, y_input+(y_input_diff+y_input_height)*11-2, x_input_width-20, y_input_height]); 
                           
%     input_fetal_4 = uicontrol(panel_fetal_params, 'Style', 'edit' ...
%                                ,'Position',[x_input, y_input+(y_input_diff+y_input_height)*10, x_input_width, y_input_height] ...
%                                , 'KeyPressFcn',@cb_enter_key); 
    popup_fetal_ftypeacc = uicontrol(panel_fetal_params,'Style','popupmenu',...
                    'String',ftypeacc_strings, ...
                    'Value',1, ...
                    'Position',[x_input, y_input+(y_input_diff+y_input_height)*10, x_input_width, y_input_height], ...
                    'Callback',@cb_popup_fetal_ftypeacc);
    label_fetal_4 = uicontrol(panel_fetal_params,'Style','text'...
                               ,'String','ftypeacc'...
                               ,'fontsize',fontSize ...    
                               ,'HorizontalAlignment','left'...              
                               ,'Position',[x_offset, y_input+(y_input_diff+y_input_height)*10-2, x_input_width-20, y_input_height]); 
                           
%     input_fetal_5 = uicontrol(panel_fetal_params, 'Style', 'edit' ...
%                                ,'Position',[x_input, y_input+(y_input_diff+y_input_height)*9, x_input_width, y_input_height] ...
%                                , 'KeyPressFcn',@cb_enter_key);
    chk_fetal_fecbt = uicontrol(panel_fetal_params, 'Style', 'checkbox' ...
                               ,'String', '' ...
                               ,'Position',[x_input+4/3*x_input_width_3+buf, y_input+(y_input_diff+y_input_height)*9, x_input_width_3, y_input_height] ...
                               , 'Callback',@cb_chk_fetal_fecbt);
    label_fetal_5 = uicontrol(panel_fetal_params,'Style','text'...
                               ,'String','fectb'...
                               ,'fontsize',fontSize ...    
                               ,'HorizontalAlignment','left'...              
                               ,'Position',[x_offset, y_input+(y_input_diff+y_input_height)*9-2, x_input_width-20, y_input_height]); 
                           
    input_fetal_6 = uicontrol(panel_fetal_params, 'Style', 'edit' ...
                               ,'Position',[x_input, y_input+(y_input_diff+y_input_height)*8, x_input_width, y_input_height] ...
                               , 'KeyPressFcn',@cb_enter_key); 
    label_fetal_6 = uicontrol(panel_fetal_params,'Style','text'...
                               ,'String','fres'...
                               ,'fontsize',fontSize ...
                               ,'HorizontalAlignment','left'...
                               ,'Position',[x_offset, y_input+(y_input_diff+y_input_height)*8-2, x_input_width-20, y_input_height]); 
                                 
%     input_fetal_7 = uicontrol(panel_fetal_params, 'Style', 'edit' ...
%                                ,'Position',[x_input, y_input+(y_input_diff+y_input_height)*7, x_input_width, y_input_height] ...
%                                , 'KeyPressFcn',@cb_enter_key); 
%     label_fetal_7 = uicontrol(panel_fetal_params,'Style','text'...
%                                ,'String','faccmean'...
%                                ,'fontsize',fontSize ...
%                                ,'HorizontalAlignment','left'...
%                                ,'Position',[x_offset, y_input+(y_input_diff+y_input_height)*7-2, x_input_width-20, y_input_height]); 
%                                  
%     input_fetal_8 = uicontrol(panel_fetal_params, 'Style', 'edit' ...
%                                ,'Position',[x_input, y_input+(y_input_diff+y_input_height)*6, x_input_width, y_input_height] ...
%                                , 'KeyPressFcn',@cb_enter_key); 
%     label_fetal_8 = uicontrol(panel_fetal_params,'Style','text'...
%                                ,'String','faccstd'...
%                                ,'fontsize',fontSize ...
%                                ,'HorizontalAlignment','left'...
%                                ,'Position',[x_offset, y_input+(y_input_diff+y_input_height)*6-2, x_input_width-20, y_input_height]); 
                                 
%     input_fetal_9 = uicontrol(panel_fetal_params, 'Style', 'edit' ...
%                                ,'Position',[x_input, y_input+(y_input_diff+y_input_height)*7, x_input_width, y_input_height] ...
%                                , 'KeyPressFcn',@cb_enter_key); 
    popup_fetal_ftraj = uicontrol(panel_fetal_params,'Style','popupmenu',...
                    'String',ftraj_strings, ...
                    'Value',1, ...
                    'Position',[x_input, y_input+(y_input_diff+y_input_height)*7, x_input_width, y_input_height], ...
                    'Callback',@cb_popup_fetal_ftraj);
    label_fetal_9 = uicontrol(panel_fetal_params,'Style','text'...
                               ,'String','ftraj'...
                               ,'fontsize',fontSize ...
                               ,'HorizontalAlignment','left'...
                               ,'Position',[x_offset, y_input+(y_input_diff+y_input_height)*7-2, x_input_width-20, y_input_height]); 
                                 
    input_fetal_10 = uicontrol(panel_fetal_params, 'Style', 'edit' ...
                               ,'Position',[x_input, y_input+(y_input_diff+y_input_height)*6, x_input_width, y_input_height] ...
                               , 'KeyPressFcn',@cb_enter_key); 
    label_fetal_10 = uicontrol(panel_fetal_params,'Style','text'...
                               ,'String','fvcg'...
                               ,'fontsize',fontSize ...
                               ,'HorizontalAlignment','left'...
                               ,'Position',[x_offset, y_input+(y_input_diff+y_input_height)*6-2, x_input_width-20, y_input_height]); 
                                 
    % List of foetuses. These can be added and removed using buttons add foetus and remove foetus. Adding a foetus will create an entry in this list.
    % When the foetus is selected, its parameters should be visible.
    list_foetus = uicontrol(panel_fetal_params,'style','list',...
                     'unit','pix',...
                     'position',[10 50 245 90],...
                     'min',0,'max',2,...
                     'fontsize',fontSize,...
                     'Callback', @cb_list_foetus, ...
                     'string',{'foetus 1', 'foetus 2'});

    % Add foetus button
    bt_add_foetus = uicontrol(panel_fetal_params...
                    ,'Style','pushbutton'...
                    ,'String','Add foetus'...
                    ,'Position',[10 10 100 30]...
                    ,'Callback',@cb_add_foetus ...
            );

    % Remove foetus button
    bt_remove_foetus = uicontrol(panel_fetal_params...
                    ,'Style','pushbutton'...
                    ,'String','Remove foetus'...
                    ,'Position',[120 10 100 30]...
                    ,'Callback',@cb_remove_foetus ...
            ); 

        
% List of scenarios for loading scenario parameters -- either run by button press, or create callback on click of list item                       
list_scenarios = uicontrol(panel_custom_scenarios,'style','list',...
                 'unit','pix',...
                 'position',[10 45 245 120],...
                 'min',0,'max',1,...
                 'fontsize',fontSize,...
                 'Callback', @cb_list_scenarios, ...
                 'string',dscen_labels);
             
    bt_import_custom = uicontrol(panel_custom_scenarios...
                    ,'Style','pushbutton'...
                    ,'String','Load'...
                    ,'Position',[10 5 100 30]...
                    ,'Callback',@cb_bt_import_custom ); 
    bt_export_custom = uicontrol(panel_custom_scenarios...
                    ,'Style','pushbutton'...
                    ,'String','Save'...
                    ,'Position',[150 5 100 30]...
                    ,'Callback',@cb_bt_export_custom); 
    
        
    % Help buttons
    bt_custom_general_help = uicontrol(panel_general_params...
                    ,'Style','pushbutton'...
                    ,'String','?'...
                    ,'Fontweight', 'bold' ...
                    ,'Position',[251 104 15 15]...
                    ,'Callback',@cb_bt_custom_general_help ...
                    ,'tooltip','Help' );
    bt_custom_fetal_help = uicontrol(panel_fetal_params...
                    ,'Style','pushbutton'...
                    ,'String','?'...
                    ,'Fontweight', 'bold' ...
                    ,'Position',[251 362 15 15]...
                    ,'Callback',@cb_bt_custom_fetal_help ...
                    ,'tooltip','Help' );
    bt_custom_noise_help = uicontrol(panel_noise_params...
                    ,'Style','pushbutton'...
                    ,'String','?'...
                    ,'Fontweight', 'bold' ...
                    ,'Position',[250 233 15 15]...
                    ,'Callback',@cb_bt_custom_noise_help ...
                    ,'tooltip','Help');
    bt_custom_mother_help = uicontrol(panel_mother_params...
                    ,'Style','pushbutton'...
                    ,'String','?'...
                    ,'Fontweight', 'bold' ...
                    ,'Position',[250 233 15 15]...
                    ,'Callback',@cb_bt_custom_mother_help ...
                    ,'tooltip','Help');
    bt_custom_controls_help = uicontrol(panel_custom_controls...
                    ,'Style','pushbutton'...
                    ,'String','?'...
                    ,'Fontweight', 'bold' ...
                    ,'Position',[251 39 15 15]...
                    ,'Callback',@cb_bt_custom_controls_help ...
                    ,'tooltip','Help' );
    bt_custom_scen_help = uicontrol(panel_custom_scenarios...
                    ,'Style','pushbutton'...
                    ,'String','?'...
                    ,'Fontweight', 'bold' ...
                    ,'Position',[251 168 15 15]...
                    ,'Callback',@cb_bt_custom_scen_help ...
                    ,'tooltip','Help');
    bt_custom_geo_help = uicontrol(panel_geometry...
                    ,'Style','pushbutton'...
                    ,'String','?'...
                    ,'Fontweight', 'bold' ...
                    ,'Position',[251 232 15 15]...
                    ,'Callback',@cb_bt_custom_geo_help ...
                    ,'tooltip','Help');        

    % Help texts    
    label_custom_general_help = uicontrol(panel_general_params_help,'Style','text'...
                               ,'String',custom_general_help_text...
                               ,'fontsize',fontSizeHelp ...        
                               ,'HorizontalAlignment','left'...  
                               ,'Visible','on' ...        
                               ,'Position',[5 5 250 105]);     
    label_custom_fetal_help = uicontrol(panel_fetal_params_help,'Style','text'...
                               ,'String',custom_fetal_help_text...
                               ,'fontsize',fontSizeHelp ...        
                               ,'HorizontalAlignment','left'...     
                               ,'Visible','on' ...        
                               ,'Position',[5 75 250 290]);     
    label_custom_noise_help = uicontrol(panel_noise_params_help,'Style','text'...
                               ,'String',custom_noise_help_text...
                               ,'fontsize',fontSizeHelp ...        
                               ,'HorizontalAlignment','left'...    
                               ,'Visible','on' ...        
                               ,'Position',[5 5 250 230]);     
    label_custom_mother_help = uicontrol(panel_mother_params_help,'Style','text'...
                               ,'String',custom_mother_help_text...
                               ,'fontsize',fontSizeHelp ...        
                               ,'HorizontalAlignment','left'...     
                               ,'Visible','on' ...        
                               ,'Position',[5 5 250 230]);     
    label_custom_controls_help = uicontrol(panel_custom_controls_help,'Style','text'...
                               ,'String',custom_controls_help_text...
                               ,'fontsize',fontSizeHelp ...        
                               ,'HorizontalAlignment','left'...     
                               ,'Visible','on' ...        
                               ,'Position',[5 5 250 40]);     
    label_custom_scen_help = uicontrol(panel_custom_scenarios_help,'Style','text'...
                               ,'String',custom_scen_help_text...
                               ,'fontsize',fontSizeHelp ...        
                               ,'HorizontalAlignment','left'...    
                               ,'Visible','on' ...        
                               ,'Position',[5 5 250 165]);     
    label_custom_geo_help = uicontrol(panel_geometry_help,'Style','text'...
                               ,'String',custom_geo_help_text...
                               ,'fontsize',fontSizeHelp ...        
                               ,'HorizontalAlignment','left'...    
                               ,'Visible','on' ...        
                               ,'Position',[5 5 250 230]);     

    % Help-Back buttons
    bt_custom_general_help_x = uicontrol(panel_general_params_help...
                    ,'Style','pushbutton'...
                    ,'String','X'...
                    ,'Fontweight', 'bold' ...
                    ,'Position',get(bt_custom_general_help, 'Position')...
                    ,'Callback',@cb_bt_custom_general_help_x ...
                    ,'tooltip','Help' ...
            );
   
    bt_custom_fetal_help_x = uicontrol(panel_fetal_params_help...
                    ,'Style','pushbutton'...
                    ,'String','X'...
                    ,'Fontweight', 'bold' ...
                    ,'Position',get(bt_custom_fetal_help, 'Position')...
                    ,'Callback',@cb_bt_custom_fetal_help_x ...
                    ,'tooltip','Help' ...
            );
      
    bt_custom_noise_help_x = uicontrol(panel_noise_params_help...
                    ,'Style','pushbutton'...
                    ,'String','X'...
                    ,'Fontweight', 'bold' ...
                    ,'Position',[250 233 15 15]...
                    ,'Callback',@cb_bt_custom_noise_help_x ...
                    ,'tooltip','Help' ...
            );
        
    bt_custom_mother_help_x = uicontrol(panel_mother_params_help...
                    ,'Style','pushbutton'...
                    ,'String','X'...
                    ,'Fontweight', 'bold' ...
                    ,'Position',[250 233 15 15]...
                    ,'Callback',@cb_bt_custom_mother_help_x ...
                    ,'tooltip','Help' ...
            );
    
    bt_custom_controls_help_x = uicontrol(panel_custom_controls_help...
                    ,'Style','pushbutton'...
                    ,'String','X'...
                    ,'Fontweight', 'bold' ...
                    ,'Position',[251 39 15 15]...
                    ,'Callback',@cb_bt_custom_controls_help_x ...
                    ,'tooltip','Help' ...
            );
    
    bt_custom_scen_help_x = uicontrol(panel_custom_scenarios_help...
                    ,'Style','pushbutton'...
                    ,'String','X'...
                    ,'Fontweight', 'bold' ...
                    ,'Position',[251 168 15 15]...
                    ,'Callback',@cb_bt_custom_scen_help_x ...
                    ,'tooltip','Help' ...
            );
    
    bt_custom_geo_help_x = uicontrol(panel_geometry_help...
                    ,'Style','pushbutton'...
                    ,'String','X'...
                    ,'Fontweight', 'bold' ...
                    ,'Position',[251 232 15 15]...
                    ,'Callback',@cb_bt_custom_geo_help_x ...
                    ,'tooltip','Help' ...
            );        
                        
% For setting enabling/disabling (custom enable, other scenarios disable)
all_custom_inputs = [input_mother_1_1;
                    input_mother_1_2;
                    input_mother_1_3;
                    input_mother_2;
                    input_mother_3;
                    popup_mother_mtypeacc;
%                     input_mother_5;
                    chk_mother_mecbt;
                    input_mother_6;
                    input_mother_7;
                    input_general_1;
                    input_general_2;
                    input_noise_1;
                    input_noise_2;
                    popup_noise_ntype;
                    input_noise_4;
                    bt_add_noise;
                    bt_remove_noise;
                    input_fetal_1_1;
                    input_fetal_1_2;
                    input_fetal_1_3;
                    input_fetal_2;
                    input_fetal_3;
                    popup_fetal_ftypeacc;
%                     input_fetal_5;
                    chk_fetal_fecbt;
                    input_fetal_6;
                    popup_fetal_ftraj
                    input_fetal_10;
                    bt_add_foetus;
                    bt_remove_foetus];        
        
        
% %% title
% 
% panel_title = uipanel('Parent', fh ...
%                     , 'Position', [0.02 0.92 0.96 0.07] ...
%                     , 'Visible', 'on' ...
%                     , 'BorderType', 'none' ...
%                     );
%                 
%       
% 
% title_label = uicontrol(panel_title,'Style','text'...
%                                ,'String',gui_title...
%                                ,'fontsize',fontSizeTitle ...        
%                                ,'HorizontalAlignment','left'...  
%                                ,'Visible','on' ...        
%                                ,'Position',[350 10 250 30]);   

%% Utility functions for MYGUI
% Updates the target_handle to display the contents of the source_handle
function new_handle = update_axes(target_ui_element, source_handle, persistent_handle)
    % target_handle is a handle to an axes element
    % source_handle is a handle to a figure element
    
    if ishandle(persistent_handle)
        all_children = num2cell(allchild(persistent_handle));
        
        % Hide all children objects followed by the parent object
        hide_elements_rec(all_children);
        set(persistent_handle,'Visible','off')
    end
    
    tmpaxes=findobj(source_handle,'Type','axes');
    new_handle = copyobj(tmpaxes,target_ui_element);
    
    % link x-axes if there is more than one plot
    axesHandles = new_handle;
    classHandles = handle(axesHandles);
    count = length(axesHandles);
    isNotInstanceOfSubtype = false(1, count);
    for i = 1:count
        isNotInstanceOfSubtype(i) = strcmp(class(classHandles(i)), 'axes') == 1;
    end
    axesHandles = axesHandles(isNotInstanceOfSubtype);
    linkaxes(axesHandles, 'x')
    
end


% recursively set children handles' 'visible' property to 'off'
function hide_elements_rec(handle)
    func = @(x) set(x, 'Visible', 'off');
    apply_to_all_children(handle, func);
end


% % recursively apply function to all children of the input handle
function apply_to_all_children(handle, func)
    if iscell(handle) && ~isempty(handle)
        for i = 1:length(handle)
            apply_to_all_children(handle{i}, func);
        end
    elseif ~iscell(handle)
         for i = 1:length(handle)
             func(handle(i));
         end
    end  
end


% Switching between the busy states
function toggle_busy(flag)
    if flag == 0
        busy_flag = 0;
        set(label_gui_busy, 'Visible', 'off');

    elseif flag == 1
        busy_flag = 1;
        set(label_gui_busy, 'Visible', 'on');
    end

end

% if flag == 1, then the makes input panels invisible and the geometry
% preview visible. also recomputes the electrode plot
% if flag == 0, it makes the input fields visible again and the preview
% plot invisible
function update_geometry_preview(flag)
   if flag == 1
        set(panel_fetal_params, 'Visible', 'off');
        set(panel_noise_params, 'Visible', 'off');
        set(panel_mother_params, 'Visible', 'off');
        set(panel_general_params, 'Visible', 'off');
        set(panel_elpos_preview, 'Visible', 'on');

        electrode_plot_handle = FECGSYN_UI_create_electrode_plot(param_struct{get(list_scenarios,'Value')});
        displayed_electrode_plot_handle = update_axes(panel_elpos_preview, electrode_plot_handle, displayed_electrode_plot_handle);
   
   elseif flag == 0
        set(panel_elpos_preview, 'Visible', 'off');
        set(panel_fetal_params, 'Visible', 'on');
        set(panel_noise_params, 'Visible', 'on');
        set(panel_mother_params, 'Visible', 'on');
        set(panel_general_params, 'Visible', 'on');
   end
end


% Adds all plots to the list on the main view and displays the first one
function show_main_plots(f, selection)
    if nargin < 1
        f = f_handles;
    end
    if nargin < 2
        selection = 1;
    end
    
    % display the first plot
    active_plot_handle = update_axes(axes_plots, f{1}.plots(1), active_plot_handle);
    
    % Update the list of plots in list_plots
    plot_titles = {};
    for i = 1:length(f)
        plot_titles{i} = f{i}.title;
    end
    plot_titles_children = {};
    for i = 1:length(f{selection}.plots)
        plot_titles_children{i} = get(f{selection}.plots(i),'name');
    end
    
    set(list_plots, 'string', plot_titles);
    set(list_plots, 'value', selection);
    set(list_plots_children, 'string', plot_titles_children);
    set(list_plots_children, 'value', 1);
    
    if length(get(list_plots_children, 'string')) > 1
        set(list_plots_children, 'Visible', 'on');
    else
        set(list_plots_children, 'Visible', 'off');
    end
end


% populates the custom view given the param struct
% Inputs: p (struct) = param as passed to run_ecg_generator
%         foetus (integer) = choice of foetus to display
%         np (integer) = choice of noise source to display
function populate_custom_view(p, foetus, ns)
    
    % Optional inputs
    if nargin < 1 || isempty(p)
        % Set foetus to 1 if there is at least one foetus
        p = param_struct{get(list_scenarios, 'Value')};
    end
    if nargin < 2 
        % Set foetus to 1 if there is at least one foetus
        if ~isempty(p.fhr);
            foetus = get(list_foetus, 'Value');
        else
            foetus = 0;
        end
    end
    if nargin < 3
        % Set ns to 1 if there is at least one noise source
        if ~isempty(p.noise_fct);
            ns = get(list_noise_sources, 'Value');
        else 
            ns = 0;
        end
    end
    
    n_foetus = length(p.fhr);
    n_ns = length(p.noise_fct);
    
    % populate foetus and noise source lists
    foetus_list = cell(n_foetus,1);
    for i = 1:n_foetus
        foetus_list{i} = sprintf('Foetus %d',i);
    end
    ns_list = cell(n_ns,1);
    for i = 1:n_ns
        ns_list{i} = sprintf('Noise Source %d',i);
    end
    set(list_foetus, 'String', foetus_list);
    set(list_noise_sources, 'String', ns_list);
    
    
    % populate general params
    set(input_general_1, 'String', num2str(p.n));
    set(input_general_2, 'String', num2str(p.fs));
    %set(input_general_3, 'String', 'elpos: [34x3 double]');
    
    
    % populate fetal params
    if foetus > 0
        % populate foetus params
        set(input_fetal_1_1, 'String', num2str(p.fheart{foetus}(1)));
        set(input_fetal_1_2, 'String', num2str(p.fheart{foetus}(2)));
        set(input_fetal_1_3, 'String', num2str(p.fheart{foetus}(3)));
        set(input_fetal_2, 'String', num2str(p.fhr(foetus)));
        set(input_fetal_3, 'String', num2str(p.facc(foetus)));
        %set(input_fetal_4, 'String', p.ftypeacc{foetus});
        idx = find(strcmp(param_struct{get(list_scenarios, 'Value')}.ftypeacc{foetus}, ftypeacc_strings));
        set(popup_fetal_ftypeacc, 'Value', idx);
%         set(input_fetal_5, 'String', p.fectb);
        set(chk_fetal_fecbt, 'Value', p.fectb);
        set(input_fetal_6, 'String', p.fres(foetus));
%         set(input_fetal_7, 'String', p.faccmean{foetus});
%         set(input_fetal_8, 'String', p.faccstd{foetus});
%         set(input_fetal_9, 'String', p.ftraj{foetus});
        idx = find(strcmp(param_struct{get(list_scenarios, 'Value')}.ftraj{foetus}, ftraj_strings));
        set(popup_fetal_ftraj, 'Value', idx);
        set(input_fetal_10, 'String', p.fvcg(foetus));
    else
        set(input_fetal_1_1, 'String', 'N/A');
        set(input_fetal_1_2, 'String', 'N/A');
        set(input_fetal_1_3, 'String', 'N/A');
        set(input_fetal_2, 'String', 'N/A');
        set(input_fetal_3, 'String', 'N/A');
        %set(input_fetal_4, 'String', 'N/A');
        set(popup_fetal_ftypeacc, 'Value', 1);
%         set(input_fetal_5, 'String', 'N/A');
        set(chk_fetal_fecbt, 'Value', 0);
        set(input_fetal_6, 'String', 'N/A');
%         set(input_fetal_7, 'String', 'N/A');
%         set(input_fetal_8, 'String', 'N/A');
%         set(input_fetal_9, 'String', 'N/A');
        set(popup_fetal_ftraj, 'Value', 1);
        set(input_fetal_10, 'String', 'N/A');
    end
    
    
    % populate noise params
    if ns > 0
        %set(input_noise_3, 'String', p.ntype{ns});
        idx = find(strcmp(param_struct{get(list_scenarios, 'Value')}.ntype{ns}, ntype_strings));
        set(popup_noise_ntype, 'Value', idx);
        set(input_noise_4, 'String', p.noise_fct_str{ns});
%         idx = find(strcmp(num2str(param_struct{get(list_scenarios, 'Value')}.noise_fct{ns}), noise_fct_strings));
%         set(popup_noise_fct, 'Value', idx);
    else
%         set(input_noise_1, 'String', 'N/A');
%         set(input_noise_2, 'String', 'N/A');
        %set(input_noise_3, 'String', 'N/A');
        set(popup_noise_ntype, 'Value', 1);
        set(input_noise_4, 'String', 'N/A');
        %set(popup_noise_fct, 'Value', 1);
    end
    set(input_noise_1, 'String', num2str(p.SNRfm));
    set(input_noise_2, 'String', num2str(p.SNRmn));
    
    % populate mother params
    set(input_mother_1_1, 'String', num2str(p.mheart(1)));
    set(input_mother_1_2, 'String', num2str(p.mheart(2)));
    set(input_mother_1_3, 'String', num2str(p.mheart(3)));
    set(input_mother_2, 'String', num2str(p.mhr));
    set(input_mother_3, 'String', num2str(p.macc));
    %set(input_mother_4, 'String', p.mtypeacc);
    idx = find(strcmp(param_struct{get(list_scenarios, 'Value')}.mtypeacc, mtypeacc_strings));
    set(popup_mother_mtypeacc, 'Value', idx);
%     set(input_mother_5, 'String', num2str(p.mectb));
    set(chk_mother_mecbt, 'Value', p.mectb);
    set(input_mother_6, 'String', num2str(p.mres));
    set(input_mother_7, 'String', num2str(p.evcg));
    
    
    % populate geometry params
    
    % hide all elements first. the required number of gui elements will be
    % unhidden depending on how many elpos valus are to be displayed
    for ii = 1:size(geo_ui_elements, 1)
        for jj = 1:size(geo_ui_elements, 2)
            set(geo_ui_elements(ii, jj), 'Visible', 'off');
        end
    end
    
    % extract the relevant elpos values
    % page min max
    %  1	1	7
    %  2	8	14
    %  3	15	28
    %  4	29	34(!)
    elpos_min = 1 + (elpos_page-1)*elpos_page_size;
    elpos_max = min([elpos_page*elpos_page_size , size(p.elpos, 1)]);
    elpos_displayed_idx = elpos_min:elpos_max;
    elpos_disp = p.elpos(elpos_min:elpos_max, :);
    
    % set visibility depending on page
    visible_electrodes = size(elpos_disp,1);
    for ii = 1:visible_electrodes
        for jj = 1:3 % the label and the three editboxes
            set(geo_ui_elements(ii,jj), 'Visible', 'on');
        end
    end
   
    % populate gui elements
    for ii = 1:size(elpos_disp, 1) % rows
        % editboxes
        % middle numerical value (r) cannot be edited
        set(geo_ui_elements(ii,2), 'String', num2str(elpos_disp(ii,1)));
        set(geo_ui_elements(ii,3), 'String', num2str(elpos_disp(ii,3)));
        
        % labels
        curr_label = ['Electrode ', num2str(elpos_displayed_idx(ii))];
        set(geo_ui_elements(ii,1), 'String', curr_label);
    end
    
    % if the preview button has been selected
    if strcmp('Hide', get(bt_geo_preview, 'String'))
        %hide_elements_rec(panel_elpos_preview);
        %electrode_plot_handle = update_axes(panel_elpos_preview, electrode_plot_handle, electrode_plot_handle);
        update_geometry_preview(0);
        update_geometry_preview(1);
    end
    
    
    drawnow();
end


% populates the custom view given the param struct
% Inputs: p (struct) = param as passed to run_ecg_generator
%         foetus (integer) = choice of foetus to display
%         np (integer) = choice of noise source to display
function save_custom_params(validate, foetus_choice, ns_choice)
% On opening the custom view, the user is presented with the values of a
% default scenario. When the user chooses to edit the scenario of choice,
% the scenario's parameter struct is copied to the last entry in
% param_struct and then the values of this struct are 
% displayed. Any user inputs will be saved to the custom param struct when
% one of the save or confirm buttons is pressed.
%
% the foetus and ns inputs are used for when the user changes the foetus or
% ns list selection and we want to save the values of the old foetus or ns.
    % The 'validate' flag checks whether the input should be checked by
    % validate_input()
    if nargin < 1
        validate = 0;
    end
    if nargin < 2
        foetus_choice = get(list_foetus, 'Value');
    end
    if nargin < 3
        ns_choice = get(list_noise_sources, 'Value');
    end
    temp = param_struct{end};

        % save general params
        temp.n = str2double(get(input_general_1, 'String'));
        temp.fs = str2double(get(input_general_2, 'String'));
        %tmp.elpos = assert_num(get(input_general_3, 'String'));
        
        % Check if any foetuses exist
        if ~isempty(get(list_foetus, 'String'))
            % save fetal params
            %temp.fheart = cell(length(get(list_foetus, 'String')), 1);
            temp.fheart{foetus_choice} = -1 * ones(1,3);
            temp.fheart{foetus_choice}(1) = str2double(get(input_fetal_1_1, 'String'));
            temp.fheart{foetus_choice}(2) = str2double(get(input_fetal_1_2, 'String'));
            temp.fheart{foetus_choice}(3) = str2double(get(input_fetal_1_3, 'String'));

            temp.fhr(foetus_choice) = str2double(get(input_fetal_2, 'String'));
            temp.facc(foetus_choice) = str2double(get(input_fetal_3, 'String'));
            %temp.ftypeacc{foetus_choice} = get(input_fetal_4, 'String');
            temp.ftypeacc{foetus_choice} = ftypeacc_strings{get(popup_fetal_ftypeacc, 'Value')};
            temp.fres(foetus_choice) = str2double(get(input_fetal_6, 'String'));
            %temp.faccmean{foetus_choice} = str2double(get(input_fetal_7, 'String'));
            %temp.faccstd{foetus_choice} = str2double(get(input_fetal_8, 'String'));
            temp.faccmean{foetus_choice} = 0;
            temp.faccstd{foetus_choice} = 1;
            temp.ftraj{foetus_choice} = ftraj_strings{get(popup_fetal_ftraj, 'Value')};
%             temp.fectb = str2double(get(input_fetal_5, 'String'));
            temp.fectb = get(chk_fetal_fecbt, 'Value');
            temp.fvcg(foetus_choice) = str2double(get(input_fetal_10, 'String')); 
            
        end
        
        
        % save noise params
        % Check if any noise sources exist
        if ~isempty(get(list_noise_sources, 'String'))
            %ns_choice = get(list_noise_sources, 'Value');
            temp.SNRfm = str2double(get(input_noise_1, 'String'));
            temp.SNRmn = str2double(get(input_noise_2, 'String'));
            %temp.ntype{ns_choice} = get(input_noise_3, 'String');
            temp.ntype{ns_choice} = ntype_strings{get(popup_noise_ntype, 'Value')};
            temp.noise_fct_str{ns_choice} = get(input_noise_4, 'String');
            temp.noise_fct{ns_choice} = eval_noise_fct(temp);  % evaluate the noise_fct string
            %temp.noise_fct{ns_choice} = noise_fct_strings{get(popup_noise_fct, 'Value')};
        end
        
        
        % save mother params
        temp.mheart = -1 * ones(1,3);
        temp.mheart(1) = str2double(get(input_mother_1_1, 'String'));
        temp.mheart(2) = str2double(get(input_mother_1_2, 'String'));
        temp.mheart(3) = str2double(get(input_mother_1_3, 'String'));
        temp.mhr = str2double(get(input_mother_2, 'String'));
        temp.macc = str2double(get(input_mother_3, 'String'));
        %temp.mtypeacc = get(input_mother_4, 'String');
        temp.mtypeacc = mtypeacc_strings{get(popup_mother_mtypeacc, 'Value')};
%         temp.mectb = str2double(get(input_mother_5, 'String'));
        temp.mectb = get(chk_mother_mecbt, 'Value');
        temp.mres = str2double(get(input_mother_6, 'String'));
        temp.evcg = str2double(get(input_mother_7, 'String'));

        
        % save geometry params
        elpos_min = 1 + (elpos_page-1)*elpos_page_size;
        elpos_max = min([elpos_page*elpos_page_size , size(temp.elpos, 1)]);
        elpos_displayed_idx = elpos_min:elpos_max;
        
        for ii = 1:length(elpos_displayed_idx)
            new_elpos_vals = zeros(3,1);
            gui_elmnts = geo_ui_elements(ii, 2:3);
            new_elpos_vals(1) = str2double(get(gui_elmnts(1), 'String'));
            new_elpos_vals(2) = temp.elpos(ii,2); % radius is always unchanged
            new_elpos_vals(3) = str2double(get(gui_elmnts(2), 'String'));
            temp.elpos(elpos_displayed_idx(ii), 1:3) = new_elpos_vals;
        end
        
        
    % Checks if all input values are in a realistic range. If so, then the
    % custom parameter structure is overwritten.
%     if validate
%         [temp, valid] = validate_param(temp);
%     else 
%         valid = 1;
%     end
    
    if 1
        param_struct{end} = temp;
    else
        disp('Error when saving custom parameters');
    end
    
end


% Checks number of foetuses and noise sources as well as the page of the
% electrodes panel so that only relevant buttons are displayed
function toggle_custom_buttons()
    % Make the add and remove foetus and noise source buttons active only if
    % the custom scenario is selected
    
    % If a default scenario is selected
    if get(list_scenarios, 'Value') < length(param_struct)
        set(bt_add_foetus, 'Visible', 'off');
        set(bt_remove_foetus, 'Visible', 'off');
        set(bt_add_noise, 'Visible', 'off');
        set(bt_remove_noise, 'Visible', 'off');
        
    elseif get(list_scenarios, 'Value') == length(param_struct)
        set(bt_add_foetus, 'Visible', 'on');
        set(bt_add_noise, 'Visible', 'on');
        
        num_foetus = length(get(list_foetus, 'String'));
        if num_foetus == 0
            set(bt_remove_foetus, 'Visible', 'off');
        else
            set(bt_remove_foetus, 'Visible', 'on');
        end
        num_ns = length(get(list_noise_sources, 'String'));
        if num_ns == 0
            set(bt_remove_noise, 'Visible', 'off');
        else
            set(bt_remove_noise, 'Visible', 'on');
        end
        
    else
        disp('Oops in toggle_add_remove_buttons()');
    end
    
    % Show the elpos next and previous buttons only if applicable
    curr_scen = get(list_scenarios, 'Value');
    curr_elpos = param_struct{curr_scen}.elpos;
    if max(elpos_displayed_idx) < length(curr_elpos) 
        % if at least one more value remains to be displayed
        set(bt_geo_next, 'Visible', 'on')
    else
        set(bt_geo_next, 'Visible', 'off')
    end
    
    if min(elpos_displayed_idx) == 1
        % if the first displayed value is the first value in the matrix
        set(bt_geo_previous, 'Visible', 'off')
    else
        set(bt_geo_previous, 'Visible', 'on')
    end
end

% takes in the param struct (in) and evaluates the noise_fct_str string at
% the time values defined by param.n
function nf = eval_noise_fct(p)
    n = p.n; % num of samples
    fs = p.fs; % frequency
    t_max = n/fs; % upper limit of time, starting at t=0
    ns = get(list_noise_sources, 'Value');
    %t = linspace(0, t_max, n);
    t = 1/fs:1/fs:(n/fs);
    nf = eval(p.noise_fct_str{ns});
    
end

%% Callbacks for MYGUI. 
% The callbacks are nested so as to allow them access to the namespace 
% of the GUI

% Callback for setting the default scenario choice via the popup menu
function cb_default_scenario_popup_menu(hObject, eventdata)
    dscen_choice = get(hObject,'Value');
    disp(['Choosing scenario: ', dscen_labels{dscen_choice}]);

end

function cb_list_plots(hObject, eventdata)
    if ~isempty(get(list_plots, 'String'))
        % Get the tag of the selected radio button. This is used to decide
        % which plot is to be displayed on the plotting axes
        curr_obj = get(list_plots,'Value');

        plot_titles_children = {};
        for i = 1:length(f_handles{curr_obj}.plots)
            plot_titles_children{i} = get(f_handles{curr_obj}.plots(i),'name');
        end
        
        set(list_plots_children, 'string', plot_titles_children);
        set(list_plots_children, 'value', 1);
        
        cb_list_plots_children([],[]);
        if length(get(list_plots_children, 'string')) > 1
            set(list_plots_children, 'Visible', 'on');
        else
            set(list_plots_children, 'Visible', 'off');
        end
    end
end

function cb_list_plots_children(hObject, eventdata)
    if ~isempty(get(list_plots_children, 'String'))
        % Get the tag of the selected radio button. This is used to decide
        % which plot is to be displayed on the plotting axes
        curr_parent = get(list_plots,'Value');
        curr_obj = get(list_plots_children,'Value');

        active_plot_handle = update_axes(axes_plots, f_handles{curr_parent}.plots(curr_obj), active_plot_handle);
    end
end

% Callback for the run button. For now it should read in the user's
% choice and run the correct scenario
function cb_run_button(hObject, eventdata)
    
    % Close any old figures in the background before generating new ones
    if ~isempty(f_handles)
        for ii = 1:length(f_handles)
            for jj = 1:length(f_handles{ii}.plots)
                close(f_handles{ii}.plots(jj))
            end
        end
    end
    
    % Remove all options in the list box to indicate that the GUI is busy
    set(list_plots, 'String', {});
    set(list_plots_children, 'String', {});
    
    % Clear the plotting axes and display the 'Busy ...' label
    if ishandle(active_plot_handle)
        all_children = num2cell(allchild(active_plot_handle));
        
        % Hide all children objects followed by the parent object
        hide_elements_rec(all_children);
        set(active_plot_handle,'Visible','off')
    end
    
    toggle_busy(1);
    
    % redraw the gui to show user that the gui is busy
    drawnow()

    if gui_mode
        % convert from popup menu index to scenario choice
        dscen_choice = get(popup_default_scenario, 'Value');
        %disp(['Running run_ecg_generator() with scenario ', num2str(dscen_choice)]);
        param = param_struct{ dscen_choice };
        out = run_ecg_generator(param,0);

%         cmqrs = adjust_mqrs_location(out.mixture(CH_CANC,:),out.mqrs,param.fs,0);
%         [res, mecg_f_handle] = mecg_cancellation(cmqrs,out.mixture(CH_CANC,:),'TS-CERUTTI',20,2,1000,debug);
%         [qrs_det,~,~] = qrs_detect(res,THR,0.150,param.fs,[],[],debug);
%         stats(out.fqrs{1}/param.fs,qrs_det/param.fs,0.05,0.5,out.param.n/param.fs,param.fs);

        % Update the figure handles
        [h1, h2, h3, h4] = FECGSYN_UI_create_fecg_plots( out, [1, 1, 1, 1], CH_CANC );
        f_handles = [h1, h2, h3, h4];
        
    else
%         % for debugging the plot display
%         f_handles(1) = figure('name','DEBUG FIGURE 1');
%         set(f_handles(1), 'Visible', 'off');
%         x = linspace(0,2*pi,200);
%         y1 = sin(x);
%         y2 = cos(x);
%         subplot(2,1,1);
%         plot(x,y1);
%         xlim([0, 2*pi]);
%         subplot(2,1,2);
%         plot(x,y2);
%         xlim([0, 2*pi]);
% 
%         f_handles(2) = figure('name','DEBUG FIGURE 2');
%         set(f_handles(2), 'Visible', 'off');
%         x = linspace(0,2*pi,200);
%         y = tan(x);
%         plot(x,y,'*');
%         xlim([0, 2*pi]);
        %ld = load('out_respiration.mat');
        %ld = load('out_noise.mat');
        %ld = load('out_noise_multiple_preg.mat');
        ld = load('out_noise_multiple_preg_small.mat');
        out = ld.out;
        [h1, h2, h3, h4] = FECGSYN_UI_create_fecg_plots( out , [1, 1, 1, 1], CH_CANC);
        f_handles = [h1, h2, h3, h4];
    end

    toggle_busy(0);
    show_main_plots(f_handles);
    
    disp('Done');
    
end

% main view import callback. imports the "out" structure
% TODO: Add param struct into the param_struct list just before 'Custom'
function cb_import_button(hObject,eventdata)
    try
        toggle_busy(1);
        [FileName,PathName,FilterIndex] = uigetfile('*.mat');
        full_path = strcat(PathName, FileName);
        imported = load(full_path);
        try 
            % This is the name of the out structure that is exported
            % with this GUI
            temp = imported.out_io;
        catch
            % Fallback name
            temp = imported.out;
        end
        out = temp;
        
        % Save the "out" param to the custom one. This may be changed in
        % the future to add an additional entry in the scenario list with
        % the same name as the file that was imported
        param_struct{end} = out.param;
        set(popup_default_scenario, 'Value', length(param_struct));
        
        [h1, h2, h3, h4] = FECGSYN_UI_create_fecg_plots( out );
        f_handles = [h1, h2, h3, h4];
        show_main_plots(f_handles);
        
        toggle_busy(0);
        msgbox('Import complete');
    catch
        disp('An error occurred while importing the results structure');
        toggle_busy(0);
    end

end


% exports the "out" structure to a mat file
function cb_export_button(hObject,eventdata)
    [FileName,PathName] = uiputfile('*.mat');
    full_path = strcat(PathName,FileName);
    if ~isempty(fieldnames(out))
        out_io = out;
        save(full_path, 'out_io');
        msgbox(['Successfully saved ', FileName]);
    else
        msgbox('You must run a scenario first!');
    end
end

% Callback for opening custom view
function cb_open_custom(hObject,eventdata)
    set(list_scenarios, 'Value', length(param_struct));
    % Copy the scenario that is to be customized to the custom slot
    param_struct{end} = param_struct{get(popup_default_scenario, 'Value')};
    
    cb_list_scenarios([],[]);
    % make sure the add/remove buttons are in the correct state
    toggle_custom_buttons();
    
    % switch to the custom view
    set(fh_main, 'Visible', 'off');
    set(fh_custom, 'Visible', 'on');
end

% Callback for confirming parameter choices
function cb_save_edit_custom(hObject,eventdata)
    curr_str = get(bt_save_edit_custom, 'String');
    
    if strcmp(curr_str, 'Edit')
        % Do the switching only is the custom scenario isn't already active
        if get(list_scenarios, 'Value') ~= length(param_struct)
            % Display chosen scenario and switch list to custom
            chosen_scen = get(list_scenarios, 'Value');
            % Copy chosen scenario to custom param structure. This will be
            % further customized by the user.
            param_struct{end} = param_struct{chosen_scen};
            % Make sure the title is retained
            param_struct{end}.title = 'Custom';

            % Make custom scenario the active one
            set(list_scenarios, 'Value', length(param_struct));
            % load the custom scenario (which is a copy of the chosen scenario)
            cb_list_scenarios(hObject, eventdata);
        
        % Switch the middle button's behaviour to 'save'
        %set(bt_save_edit_custom, 'String', 'Save');
        end
        
    elseif strcmp(curr_str, 'Save')
        % save all inputs to the custom scenario and validate input
        save_custom_params(1, selected_foetus, selected_ns);
        
    else
        disp('Oops in cb_save_edit_custom()');
    end
end

% Callback for confirming parameter choices
function cb_confirm_custom(hObject,eventdata)
   save_custom_params();
   
   set(fh_main, 'Visible', 'on');
   set(fh_custom, 'Visible', 'off');
   
   % Set the selected value of the popup menu to 'Custom'
   set(popup_default_scenario, 'Value', length(param_struct));
end

% Callback for cancelling custom inputs
function cb_back_custom(hObject,eventdata)
    % If the custom scenario is active, then save the latest parameter
    % values 
%     if get(list_scenarios, 'Value') == length(param_struct)
%         save_custom_params(0, selected_foetus, selected_ns);
%     end
    
    % Set the selected value of the popup menu to the selected scenario
    set(popup_default_scenario, 'Value', get(list_scenarios,'Value'));
    
    % change face back to front
    set(fh_main, 'Visible', 'on');
    set(fh_custom, 'Visible', 'off');
end

% Callback for adding foetus
function cb_add_foetus(hObject,eventdata)
    % Save custom parameters inputted so far
    save_custom_params(0, selected_foetus, selected_ns);
    
    % Increase number of foetuses in the custom param structure
    foetus_str = get(list_foetus, 'String');
    old_num_foetus = length(foetus_str);
    new_num_foetus = old_num_foetus + 1;
    foetus_str{new_num_foetus} = sprintf('foetus %d', new_num_foetus);
    set(list_foetus, 'String', foetus_str);
    set(list_foetus, 'Value', new_num_foetus);
    
    set(input_fetal_1_1, 'String', -0.31416);
    set(input_fetal_1_2, 'String', 0.35);
    set(input_fetal_1_3, 'String', -0.3);
    set(input_fetal_2, 'String', 150);
    set(input_fetal_3, 'String', 0);
    %set(input_fetal_4, 'String', 'none');
    set(popup_fetal_ftypeacc, 'Value', 1);
    %set(input_fetal_5, 'String', 0);
    set(chk_fetal_fecbt, 'Value', 0);
    set(input_fetal_6, 'String', 0);
%     set(input_fetal_7, 'String', 0);
%     set(input_fetal_8, 'String', 1);
%     set(input_fetal_9, 'String', 'none');
    set(popup_fetal_ftraj, 'Value', 1);
    set(input_fetal_10, 'String', fVCG);
    
    selected_foetus = new_num_foetus;
    
    % Save new foetus in case user switches foetus immediately
    save_custom_params(0, selected_foetus, selected_ns);
    
    % Make sure the buttons are displayed correctly after this change
    toggle_custom_buttons();

end

% Callback for removing foetus
function cb_remove_foetus(hObject,eventdata)
    foetus_choice = get(list_foetus, 'Value');

    temp = param_struct{end};
    
    temp.fheart(foetus_choice) = [];
    temp.fhr(foetus_choice) = [];
    temp.facc(foetus_choice) = [];
    temp.ftypeacc{foetus_choice} = [];
    temp.fres(foetus_choice) = [];
    %temp.fectb = [];
    
    % Update custom param struct
    param_struct{end} = temp;
    
    % Update the list of foetuses
    if foetus_choice > 1
        selected_foetus = foetus_choice - 1;
    else
        selected_foetus = 1;
    end
    foetus_str = get(list_foetus, 'String');
    foetus_str(end) = [];
    set(list_foetus, 'Value', selected_foetus);
    set(list_foetus, 'String', foetus_str);
    populate_custom_view();
    
    % Make sure the buttons are displayed correctly after this change
    toggle_custom_buttons();
end

% Callback for adding noise
function cb_add_noise(hObject,eventdata)
    % Save custom parameters inputted so far
    save_custom_params(0, selected_foetus, selected_ns);
    
    %set(input_noise_3, 'String', 'MA');
    set(popup_noise_ntype, 'Value', 1);
    set(input_noise_4, 'String', 1);
    %set(popup_noise_fct, 'Value', 1);
        
    % Increase number of noise in the custom param structure
    noise_str = get(list_noise_sources, 'String');
    old_num_noise = length(noise_str);
    new_num_noise = old_num_noise + 1;
    noise_str{new_num_noise} = sprintf('Noise Source %d', new_num_noise);
    set(list_noise_sources, 'String', noise_str);
    set(list_noise_sources, 'Value', new_num_noise);
    
    selected_ns = new_num_noise;
    
    % Save new foetus in case user switches foetus immediately
    save_custom_params(0, selected_foetus, selected_ns);
    
    % Make sure the buttons are displayed correctly after this change
    toggle_custom_buttons();

end

% Callback for removing noise
function cb_remove_noise(hObject,eventdata)
    
    ns_choice = get(list_noise_sources, 'Value');
    
    temp = param_struct{end};
    
    % Delete relevant data
    temp.ntype(ns_choice) = [];
    temp.noise_fct_str(ns_choice) = '';
    temp.noise_fct(ns_choice) = [];
       
    % Update custom param struct
    param_struct{end} = temp;
    
    % Update the list of foetuses
    if ns_choice > 1
        selected_ns = ns_choice - 1;
    else
        selected_ns = 1;
    end
    ns_str = get(list_noise_sources, 'String');
    ns_str(end) = [];
    set(list_noise_sources, 'Value', selected_ns);
    set(list_noise_sources, 'String', ns_str);
    
    populate_custom_view();
    
    % Make sure the buttons are displayed correctly after this change
    toggle_custom_buttons();
end

function cb_list_foetus(hObject,eventdata)
    save_custom_params( 0, selected_foetus, selected_ns);
    selected_foetus = get(list_foetus, 'Value');
    populate_custom_view([], selected_foetus, selected_ns);
end

function cb_list_noise_sources(hObject,eventdata)
    save_custom_params( 0, selected_foetus, selected_ns);
    selected_ns = get(list_noise_sources, 'Value');
    populate_custom_view([], selected_foetus, selected_ns);
end

function cb_list_scenarios(hObject,eventdata)
    
    % If the custom scenario is active, then change the middle button's
    % behaviour to save, otherwise edit
    scen_id = get(list_scenarios,'Value');
%     if scen_id == length(param_struct)
%         set(bt_save_edit_custom, 'String', 'Save');
%         %save_custom_params();
%     else
%         set(bt_save_edit_custom, 'String', 'Edit');
%     end
    
    % get the relevant param struct
    param = param_struct{scen_id};
    
    selected_foetus = ~isempty(get(list_foetus,'String'));
    selected_ns = ~isempty(get(list_noise_sources,'String'));
    elpos_page = 1;
    
    % Update add/remove buttons
    toggle_custom_buttons();
    
    if scen_id == length(param_struct)
        set(panel_geometry_no_edit, 'Visible', 'off')
        set(panel_geometry, 'Visible', 'on')
    else
        set(panel_geometry_no_edit, 'Visible', 'on')
        set(panel_geometry, 'Visible', 'off')
    end
    
    if get(list_scenarios,'value') == length(param_struct)
        for i = 1:length(all_custom_inputs)
        set(all_custom_inputs(i),'enable','on')
        end
    else
        for i = 1:length(all_custom_inputs)
        set(all_custom_inputs(i),'enable','off')
        end
    end
    
    % display the data
    populate_custom_view(param);
end

function cb_bt_geo_next(hObject,eventdata)
    save_custom_params();
    elpos_page = elpos_page + 1;
    pause(0.2)
    populate_custom_view();
    toggle_custom_buttons();
end

function cb_bt_geo_previous(hObject,eventdata)
    save_custom_params();
    pause(0.2)
    elpos_page = elpos_page - 1;
    populate_custom_view();
    toggle_custom_buttons();
end

function cb_bt_geo_preview(hObject,eventdata)
%     update_geometry_preview(0);
    update_geometry_preview(1);
    set(bt_geo_hide, 'Visible', 'on');
    set(bt_geo_hide_no_edit, 'Visible', 'on');
    set(bt_geo_refresh, 'Visible', 'on');
    set(bt_geo_preview, 'Visible', 'off');
end

function cb_bt_geo_hide(hObject,eventdata)
    set(bt_geo_hide, 'Visible', 'off');
    set(bt_geo_hide_no_edit, 'Visible', 'off');
    set(bt_geo_refresh, 'Visible', 'off');
    set(bt_geo_preview, 'Visible', 'on');
    update_geometry_preview(0);
end


function cb_bt_geo_refresh(hObject,eventdata)
    save_custom_params();
    populate_custom_view();
    update_geometry_preview(0);
    update_geometry_preview(1);
end


function cb_bt_import_custom(hObject,eventdata)
    try
        [FileName,PathName,FilterIndex] = uigetfile('*.mat');
        full_path = strcat(PathName, FileName);
        imported = load(full_path);
        try 
            % This is the name of the param structure that is exported
            % with this GUI
            temp = imported.param_io;
        catch
            % Fallback name
            temp = imported.param;
        end

        temp.title = 'Custom'; 
        param_struct{end} = temp;

    catch
        disp('An error occurred while importing parameters');
    end

   
    % populate foetus and noise source lists
    p = temp;
    if ~isempty(p.fhr);
        n_foetus = length(p.fhr);
        foetus_list = cell(n_foetus,1);
        for i = 1:n_foetus
            foetus_list{i} = sprintf('foetus %d',i);
        end
        set(list_foetus, 'String', foetus_list);
        set(list_foetus, 'Value', 1);
    else 
        foetus_list = {};
        set(list_foetus, 'String', foetus_list);
    end
    if ~isempty(p.noise_fct);
        n_ns = length(p.noise_fct);
        ns_list = cell(n_ns,1);
        for i = 1:n_ns
            ns_list{i} = sprintf('Noise Source %d',i);
        end
        set(list_noise_sources, 'String', ns_list);
        set(list_noise_sources, 'Value', 1);
    else
        ns_list = {};
        set(list_noise_sources, 'String', ns_list);
    end
    
    
    
    
    
    set(list_scenarios, 'Value', length(param_struct));
    cb_list_scenarios(hObject,eventdata);
    %disp('TODO: cb_bt_import_custom()');
end

function cb_bt_export_custom(hObject,eventdata)
    save_custom_params();
    chosen_scen = get(list_scenarios, 'Value');
    param_io = param_struct{chosen_scen};
    [FileName,PathName] = uiputfile('*.mat');
    full_path = strcat(PathName,FileName);
    save(full_path, 'param_io');
end

% Saves data, then runs the callbacks of the back button and then the
% run button on the main face of the GUI
function cb_bt_run_custom(hObject,eventdata)
    save_custom_params();
    cb_back_custom(hObject,eventdata);
    % Adding an artifial pause so that the whole of the main face can
    % be rendered before the GUI becomes less responsive due to
    % run_ecg_generator()
    drawnow(); pause(0.5);

    cb_run_button(hObject,eventdata);
    %disp('TODO: cb_bt_run_custom()');
end

function cb_enter_key(hObject,eventdata)
%     disp(['eventdata key pressed: ', eventdata.Key])
% 
    if strcmp('on', get(fh_custom, 'Visible'))
        if strcmp(eventdata.Key, 'return')
            % if the custom view is active and the key is 'return', then
            % save parameters and repopulate the gui
            save_custom_params();
%             populate_custom_view();
%             toggle_custom_buttons();
        end
    end

end

% Help button callbacks
    function cb_bt_run_help(hObject,eventdata)
        set(panel_defaults_help,'Visible','on')
        set(panel_defaults,'Visible','off')
        set(list_plots_children,'Visible','off')
    end
    function cb_bt_custom_general_help(hObject,eventdata)
        set(panel_general_params_help,'Visible','on')
        set(panel_general_params,'Visible','off')
    end

    function cb_bt_custom_fetal_help(hObject,eventdata)
        set(panel_fetal_params_help,'Visible','on')
        set(panel_fetal_params,'Visible','off')
    end

    function cb_bt_custom_noise_help(hObject,eventdata)
        set(panel_noise_params_help,'Visible','on')
        set(panel_noise_params,'Visible','off')
    end

    function cb_bt_custom_mother_help(hObject,eventdata)
        set(panel_mother_params_help,'Visible','on')
        set(panel_mother_params,'Visible','off')
    end

    function cb_bt_custom_controls_help(hObject,eventdata)
        set(panel_custom_controls_help,'Visible','on')
        set(panel_custom_controls,'Visible','off')
    end

    function cb_bt_custom_scen_help(hObject,eventdata)
        set(panel_custom_scenarios_help,'Visible','on')
        set(panel_custom_scenarios,'Visible','off')
    end

    function cb_bt_custom_geo_help(hObject,eventdata)
        set(panel_geometry_help,'Visible','on')
        set(panel_geometry,'Visible','off')
    end

% Help/Back button callbacks
    function cb_bt_run_help_x(hObject,eventdata)
        set(panel_defaults_help,'Visible','off')
        set(panel_defaults,'Visible','on')
        set(list_plots_children,'Visible','on')
    end
    function cb_bt_custom_general_help_x(hObject,eventdata)
        set(panel_general_params_help,'Visible','off')
        set(panel_general_params,'Visible','on')
    end

    function cb_bt_custom_fetal_help_x(hObject,eventdata)
        set(panel_fetal_params_help,'Visible','off')
        set(panel_fetal_params,'Visible','on')
    end

    function cb_bt_custom_noise_help_x(hObject,eventdata)
        set(panel_noise_params_help,'Visible','off')
        set(panel_noise_params,'Visible','on')
    end

    function cb_bt_custom_mother_help_x(hObject,eventdata)
        set(panel_mother_params_help,'Visible','off')
        set(panel_mother_params,'Visible','on')
    end

    function cb_bt_custom_controls_help_x(hObject,eventdata)
        set(panel_custom_controls_help,'Visible','off')
        set(panel_custom_controls,'Visible','on')
    end

    function cb_bt_custom_scen_help_x(hObject,eventdata)
        set(panel_custom_scenarios_help,'Visible','off')
        set(panel_custom_scenarios,'Visible','on')
    end

    function cb_bt_custom_geo_help_x(hObject,eventdata)
        set(panel_geometry_help,'Visible','off')
        set(panel_geometry,'Visible','on')
    end

    function cb_bt_about(hObject,eventdata)
%         msgbox(about_text, 'About');
        set(axes_plots,'Visible','off')
        set(panel_about,'Visible','on')
    end

    function cb_bt_about_close(hObject,eventdata)
%         msgbox(about_text, 'About');
        set(axes_plots,'Visible','on')
        set(panel_about,'Visible','off')
    end

    function cb_popup_mother_mtypeacc(hObject,eventdata)
        save_custom_params();
    end

    function cb_popup_noise_ntype(hObject,eventdata)
        save_custom_params();
    end

%     function cb_popup_noise_fct(hObject,eventdata)
%         save_custom_params();
%     end

    function cb_popup_fetal_ftypeacc(hObject,eventdata)
        save_custom_params();
    end

    function cb_popup_fetal_ftraj(hObject,eventdata)
        save_custom_params();
    end

    function cb_chk_fetal_fecbt(hObject,eventdata)
        save_custom_params();
    end

    function cb_popup_geo_num(hObject,eventdata)
        sel_num = geo_num_labels(get(popup_geo_num, 'Value'));
        p = param_struct{end};
        default_elpos = param_struct{1}.elpos;
        
        old_length = length(p.elpos);
        new_length = sel_num;
        
        if old_length > new_length
            % just shortening the elpos matrix
            p.elpos = p.elpos(1:new_length, :);
            param_struct{end} = p;
            
        elseif old_length < new_length
            % adding default values as padding to grow to desired size
            extra_elpos_idx = (old_length+1):new_length;
            extra_elpos = default_elpos(extra_elpos_idx, :);
            p.elpos(extra_elpos_idx, :) = extra_elpos;
            param_struct{end} = p;
        end
        
        toggle_custom_buttons()
        populate_custom_view();
        
    end



%% Make the GUI visible when everything is loaded

% display the param struct of the first scenario in the custom view
populate_custom_view(param_struct{1});

set(fh, 'Visible', 'on');


end
