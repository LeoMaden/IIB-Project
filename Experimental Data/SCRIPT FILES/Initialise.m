% Load variables and initialise measurement hardware for the red tunnel EDF rig
% 
% August 2022


%% Initialise NI 2090 using PXIe scripts for voltage readings

% Use NI PXIe 1073 for voltage readings
pxie.sesh = 'ni'; pxie.dev = [repmat({'Dev1'},[1 16])];
pxie.type = 'V'; pxie.coupling = repmat({'DC'},[1 16]); 
pxie.ai = [0:15]; pxie.nchan = 16; pxie.rate = 80e3;
pxie.time = 2; pxie.chan = 1:3;

s.pxie = exp_pxie_initialise(pxie);
disp('NI unit initialised')

%% Set default DSA logging variables
calz = 1;
dsa.ip = '129.169.103.137';
dsa.port = 23; dsa.nchan = 16;
dsa.period = 125; dsa.avg = 50;
s.dsa = exp_dsa_initialise(dsa);
disp('DSA initialised')

%% Vaisala PTH sensor
s.ptu.loc = 'HSE';
[Pa,Ta,Ha] = getCurrentPTH(s.ptu.loc);
rho = Pa/(287*Ta);
disp('PTU initialised')

%% Initialise TC08

tc08.library = {'C:\Users\lm859\Apps\UsbTC08sdk_r10_5_0_28\x64\usbtc08.dll' ...
    'C:\Users\lm859\Apps\UsbTC08sdk_r10_5_0_28\usbtc08.h'};
tc08.type = 'K'; tc08.units = 0; tc08.period = 100; tc08.nchan = 1;  tc08.chan = 1; 
s.tc08 = exp_tc08_initialise(tc08);
disp('TC08 initialised')
%% Instrument channel numbers

% Record total number of channels for all instruments


% Thermocouple Channels
N.T_motor = 1;

% Voltage channels
N.V_shaft = 1;
N.V_long = 2;
N.V_lat = 3;

% Pressure Channels
N.Po_tunnel = 1;
N.Ps_tunnel = 2;
%% Initialise Arduino Connection

%Create Arduino Object
s.arduino = arduino();

%Create Servo Object
s.servo = servo(s.arduino, "D2", 'MinPulseDuration', 1000e-6, 'MaxPulseDuration', 2000e-6);
writePosition(s.servo,0);
disp('Speed control ardunio initialised')

%% Initialse CJC stepper controller
S.COM ="COM9";
S.Vmin = [20000.,10000.,5000.]; %minimum speed step/s
S.Vmax = [20000.,15000.,10000.]; %maximum speed step/s
S.Acc = [1500.,15000.,1500.]; %acceleration step/s^2
S.Conv = [1600,1600,1600]; %conversion step per mm/deg
S.HomePos = [0.,0.,0.];%position of home switch
S.Nchannel = 3; % number of channels to use

[s.step,S]=exp_stepper_init(S);

% Traverse Channel Numbers
N.Y = 1;
N.Z = 0;

disp('Stepper Controller initialised')
%% Other Properties
% Gas constants
exp_air

