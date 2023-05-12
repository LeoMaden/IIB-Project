%% Matlab code to run red tunnel tests of EDFs


%% Initialisation


if exist('s','var') == 0
    clear; close all; clc;
    Initialise

else
    clear e
    clc
end
%% Set Jobs

calib_speed = 0;        % Calicrate EDF speed control
calib_load = 0;         % Calibrate Load Cells
run_loads = 1;          % Run loadcell experiment
run_pressures = 0;      % Run Pressure experiment
run_traverse = 0;       % Traverse Experiement
pressure_rise = 0;      % run pressure rise measurement





%% Calibration
if run_loads == 1
    if isfile('load_calib.mat') == 0 || calib_load == 1
        calibrate_loadcell
    else
        load('load_calib.mat')

        %Re-measure intercept for each new run
         %[V1,tdump] = exp_pxie_read(s.pxie);
         %load_calib(1,2) = - load_calib(1,1) * mean(V1(:,N.V_long));
         %load_calib(2,2) = - load_calib(2,1) * mean(V1(:,N.V_lat));
    end
end


%% Experiments

%e is experimental data structure

% set the geometry to run
% IMP5-V0
% e.geom.name = 'IPM5-75';
% e.geom.rhub = 0.0220;
% e.geom.rtip = 0.0718;

% % Lift Fan
e.geom.name = 'Lift_Fan';
e.geom.rhub = 0.302*0.106;
e.geom.rtip = 0.106;
% 
% % IPM5-V1
% e.geom.name = 'IPM5_V1';
% e.geom.rhub = 28.57/1000;
% e.geom.rtip = 0.08;

e.geom.rmean = sqrt((e.geom.rhub^2+e.geom.rtip^2)/2);
e.geom.A1 = pi * (e.geom.rtip^2 - e.geom.rhub^2);
%e.nose = 'IPM5-001';                                    % Nose Name for different Pressure Tapping Inserts
e.atm = struct('Pa',Pa,'Ta',Ta,'Ha',Ha);                % Add Atmospheric Data to structure
e.raw.V_supply = 44;                                    % Power Supply Voltage
e.rho = rho;
e.env = struct('Pa',Pa,'Ta',Ta,'Ha',Ha);


if run_loads == 1
    useSpeedMap = 0;
    
    e.aoa = [45 90];                                          % Set Angle of Attack range
                        
    % DO NOT GO ABOVE 0.3 FOR THE UPPER SPEED for IPM5_V0 on 44Volts
    e.tunnel_switch = [0 1 3 4 5 8 9];
    e.pwm_speed = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6];                % Set Arduino PWM range (0 is off, 1 is full speed)
    e.calib.load = load_calib;                          % Add Calibration Data
    if calib_speed == 1
        mapspeed
    end
    thrust_experiment
    close all
   
    %thrust_process
    savefile
end

if run_pressures == 1
%     addpath tapping_positions
%     
%     load('IPM5_001.mat')
%     e.coords = plot_coords;
%     e.aoa = 0;
%     e.tunnel_switch = [0 1 3 5 8 10 11];
%     e.pwm_speed = [0 0.0900 0.1017 0.1133 0.1250 0.1367 0.1483 0.1600 0.1717 0.1833 0.1950];
%     e.configs = ceil(length(e.coords(:,1))/11);
%     pressure_experiment
    e.tunnel_switch = [0 1 5 10];
    e.pwm_speed = [0 linspace(0.1,1,5)];
    Lift_Fan_Pressure_Experiment
end

if pressure_rise == 1
    e.tunnel_switch = [0 1 3 5 8 10 11];
    e.pwm_speed = [0 linspace(0.11,0.25,10)];
    rtr_pressure_rise
end

if run_traverse == 1
    e.aoa = 0;
    e.tunnel_switch = 1;
    e.pwm_speed = 0.2;
    Traverse_RTR
end

answer_dis = questdlg('Disconnect Instruments?', 'Instrument Disconnect',...
    'Yes', 'No', 'No');
% Handle response
switch answer_dis
    case 'Yes'
        Exp_Disconnect
    case 'No'        
end