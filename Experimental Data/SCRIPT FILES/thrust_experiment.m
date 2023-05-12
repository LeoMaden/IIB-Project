%% Carry Out Thrust and Drag Measurements

%% Set Up data and devices
close all;

% Pre allocate memory for speed - Experimental Data
e.raw.V_shaft = zeros(length(e.aoa), length(e.tunnel_switch), length(e.pwm_speed), (pxie.rate * pxie.time));    % Shaft Speed Signal
e.raw.t = zeros(length(e.aoa), length(e.tunnel_switch), length(e.pwm_speed), (pxie.rate * pxie.time));          % Voltage time readings
e.raw.V_long = zeros(length(e.aoa), length(e.tunnel_switch), length(e.pwm_speed), (pxie.rate * pxie.time));     % longitudinal Load zeros Signal
e.raw.V_lat = zeros(length(e.aoa), length(e.tunnel_switch), length(e.pwm_speed), (pxie.rate * pxie.time));      % Lateral Load zeros Signal
e.raw.V = zeros(length(e.aoa), length(e.tunnel_switch), length(e.pwm_speed), 3, (pxie.rate * pxie.time));       % Voltages
e.raw.pxie_time = zeros(length(e.aoa), length(e.tunnel_switch), length(e.pwm_speed), (pxie.rate * pxie.time));  % PXIe time measurements
e.raw.V_Long_av = zeros(length(e.aoa), length(e.tunnel_switch), length(e.pwm_speed));                           % Mean of Longitudinal Signal
e.raw.V_Lat_av = zeros(length(e.aoa), length(e.tunnel_switch), length(e.pwm_speed));                            % Mean of Lateral Signal
e.raw.P = zeros(length(e.aoa), length(e.tunnel_switch), length(e.pwm_speed), 16);                               % Pressures
e.raw.T_motor = zeros(length(e.aoa), length(e.tunnel_switch), length(e.pwm_speed));                             % Motor Temperature
e.raw.I_supply = zeros(length(e.aoa), length(e.tunnel_switch), length(e.pwm_speed));                            % Power Supply Current
e.raw.shaft_freq = zeros(length(e.aoa), length(e.tunnel_switch), length(e.pwm_speed));                          % Motor Speed

% Pre-allocate memory for processed results

e.long_force = zeros(length(e.aoa), length(e.tunnel_switch), length(e.pwm_speed));                              % longitudinal Force Sensor
e.lat_force = zeros(length(e.aoa), length(e.tunnel_switch), length(e.pwm_speed));                               % Lateral Force Sensor
e.J = zeros(length(e.aoa), length(e.tunnel_switch), length(e.pwm_speed));                                       % Advance Ratio

% Set Motor to be ready to run;
writePosition(s.servo, 0);
uiwait(msgbox('power on esc'));

uiwait(msgbox('Click okay when ready to run'));

%% Initialise Live Plotting

% Monitoring Window
rig_monitor = figure();
t_delay = 200; % Set scrolling time delay
% Window Layout
pos.T = [0.08 0.71 0.37 0.25];
pos.P = [0.08 0.38 0.37 0.25];
pos.Power = [0.08 0.06 0.37 0.25];
pos.Long = [0.55 0.06 0.37 0.25];
pos.Lat = [0.55 0.38 0.37 0.25];
pos.Shaft = [0.55 0.71 0.37 0.25];


a.T = axes('Position', pos.T); hold on; grid on; box on;
ylabel('Temperature (C)')
%axis([-t_delay 0 0 100]);

a.P = axes('Position', pos.P); hold on; grid on; box on;
ylabel('Pressure (Pa)')
%axis([-t_delay 0 0 10000]);

a.Long = axes('Position', pos.Long); hold on; grid on; box on;
ylabel('Longitudinal Voltage)')
%axis([-t_delay 0 -2 2]);

a.Lat = axes('Position', pos.Lat); hold on; grid on; box on;
ylabel('Lateral Voltage')
%axis([-t_delay 0 -2 2]);

a.Shaft = axes('Position', pos.Shaft); hold on; grid on; box on;
ylabel('Shaft Speed (rad/s)')
%axis([-t_delay 0 0 1000]);

a.Power = axes('Position', pos.Power); hold on; grid on; box on;
ylabel('Power (W)')


l.T = plot(a.T, nan, nan);
l.Po = plot(a.P, nan, nan);
l.Ps = plot(a.P, nan, nan);
l.Long = plot(a.Long, nan, nan);
l.Lat = plot(a.Lat, nan, nan);
l.Shaft = plot(a.Shaft, nan, nan);
l.Power = plot(a.Power,nan, nan);

hold off;

%% Test Load Cells to check calibration

% Measure zero laod
%% Record Data
for x = 1:length(e.aoa)
    uiwait(msgbox(['Set Angle of Attack to ' num2str(e.aoa(x))]))
    
    for y = 1:length(e.tunnel_switch)
        uiwait(msgbox(['Set Tunnel Speed Switch to ' num2str(e.tunnel_switch(y))]))

        check_state = figure();
        ax.temp = axes('Position', [0.08 0.56 0.84 0.38]); hold on; grid on; box on;
        ylabel('Temperature (C)')
        ax.pres = axes('Position', [0.08 0.06 0.84 0.38]); hold on; grid on; box on;
        xlabel('time (s)'); ylabel('Pressure (Pa)');
        ButtonHandle = uicontrol('Style', 'PushButton', 'String', 'Press Button Once Pressures settled', 'Callback', 'g=g+1;');
        l.temp = plot(ax.temp, nan, nan);
        l.pres0 = plot(ax.pres, nan, nan);
        l.presS = plot(ax.pres, nan, nan);
        i = 1;
        g = 0;
        pres = zeros(1,16);
        temp = nan;
        
        linkaxes([ax.temp,ax.pres],'x');
        
        time = clock;
        time_plot = 0;
        while g == 0
            pres = [pres; exp_dsa_read(s.dsa)];
            temp = [temp; exp_tc08_read(s.tc08)];
            time = [time; clock];
            time_plot = [time_plot ; etime(time(end,:),time(1,:))];
            
            set(l.temp,'xdata',time_plot,'ydata',temp);
            set(l.pres0,'xdata',time_plot,'ydata',pres(:,N.Po_tunnel));
            set(l.presS,'xdata',time_plot,'ydata',pres(:,N.Ps_tunnel));
            axis([time_plot(end)-100 time_plot(end)+100 min(pres(1:2))-10 max(pres(1:2))+100]);
            drawnow;
            pause(0.2)
        end
        clear temp pres i;
        close(check_state);
        figure(rig_monitor);
        
   
        for z = 1:length(e.pwm_speed)
            disp(['Writing Speed ' num2str(e.pwm_speed(z))])
           
            rampSpeed(s.servo, e.pwm_speed(z));                         % Set Motor Speed
            pause(2)
            
            [q, t] = exp_pxie_read(s.pxie);                                 % Read PXIe
            e.raw.t(x,y,z,:) = t;
            e.raw.V_shaft(x,y,z,:) = q(:,N.V_shaft);                  
            e.raw.V_long(x,y,z,:) = q(:,N.V_long);
            e.raw.V_lat(x,y,z,:) = q(:,N.V_lat);
            e.raw.pxie_time(x,y,z,:) = t;
            e.raw.V_Long_av(x,y,z) = mean(q(:,N.V_long));
            e.raw.V_Lat_av(x,y,z) = mean(q(:,N.V_lat));
            if useSpeedMap ~= 1
                e.raw.shaft_freq(x,y,z) = exp_shaft_freq_fft(s.pxie, 1, q(:,N.V_shaft), t);
            end
            pause(0.05)
            e.raw.T_motor(x,y,z) = exp_tc08_read(s.tc08);

            e.raw.P(x,y,z,:) = exp_dsa_read(s.dsa);                         % Read DSA
            
            %e.raw.T_motor(x,y,z) = exp_tc08_read(s.tc08);                   % Read TC08
            disp(['motor temp = ' num2str(e.raw.T_motor(x,y,z))])
%             if e.raw.T_motor(x,y,z) > 65
%                 disp('Cooling Motor')
%                 rampSpeed(s.servo,0.09)
%                 pause(60)
%             end
            
            e.raw.I_supply(x,y,z) = input('Enter Power Supply Current');    % Read Supply Current
            % enter currents at end
            % Live Plotting
            set(l.T,'xdata',1:z,'ydata',e.raw.T_motor(1:z));
            set(l.Po,'xdata',1:z,'ydata',e.raw.P(x,y,1:z,N.Po_tunnel));
            set(l.Ps,'xdata',1:z,'ydata',e.raw.P(x,y,1:z,N.Ps_tunnel));
            set(l.Long,'xdata',1:z,'ydata',e.raw.V_Long_av(1:z));
            set(l.Lat,'xdata',1:z,'ydata',e.raw.V_Lat_av(1:z));
            set(l.Shaft,'xdata',1:z,'ydata',e.raw.shaft_freq(1:z));
            set(l.Power,'xdata',1:z,'ydata',e.raw.V_supply.*e.raw.I_supply(x,y,1:z))
        end
        if strcmp(e.geom.name,'IPM5-75') == 1
            rampSpeed(s.servo,0.12)
        else
            rampSpeed(s.servo, 0)
        end
    end
end

% enter currents at end
% for x = 1:length(e.aoa)
%     for y = 1:length(e.tunnel_switch)
%         for z = 1:length(e.pwm_speed)
%             e.raw.I_supply(x,y,z) = input(['Enter Power Supply Current ' num2str(z)]);    % Read Supply Current
%         end
%     end
% end

if useSpeedMap == 1
    for i = 1:length(e.tunnel_switch)
        e.raw.shaft_freq(:,i,:) = speedmap.w;
    end
end

save('thrust_data.mat','e','N')
