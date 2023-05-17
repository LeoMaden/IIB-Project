% Sampling rate for calculating shaft frequency
f_sample = 80000;

filename = "run2";

data = load(filename + ".mat");
e = data.e;
N = data.N;

%% Fix run1
if filename == "run1"
    e.aoa(end-1:end) = [];
    e.tunnel_switch(end) = [];

    e.raw.V_long(2:3, :, :, :) = [];
    e.raw.V_long(:, end, :, :) = [];

    e.raw.V_lat(2:3, :, :, :) = [];
    e.raw.V_lat(:, end, :, :) = [];

    e.raw.P(2:3, :, :, :) = [];
    e.raw.P(:, end, :, :) = [];
    e.raw.P(:, :, :, 3:end) = [];

    e.raw.T_motor(2:3, :, :) = [];
    e.raw.T_motor(:, end, :) = [];

    e.raw.I_supply(2:3, :, :) = [];
    e.raw.I_supply(:, end, :) = [];

    e.raw.shaft_freq(2:3, :, :) = [];
    e.raw.shaft_freq(:, end, :) = [];
end

%% Extract data from e
geom = e.geom;
atm = e.atm;

V_supply = e.raw.V_supply;
I_supply = e.raw.I_supply;
P_total = e.raw.P(:, :, :, N.Po_tunnel);
P_static = e.raw.P(:, :, :, N.Ps_tunnel);
T_motor = e.raw.T_motor;

aoa = e.aoa;
tunnel_switch = e.tunnel_switch;
pwm_speed = e.pwm_speed;

% Calculate shaft frequency
f_shaft = zeros(length(aoa), length(tunnel_switch), length(pwm_speed));

for i = 1:length(aoa)
    for j = 1:length(tunnel_switch)
        for k = 1:length(pwm_speed)
            V_shaft = squeeze(e.raw.V_shaft(i, j, k, :));
            t = squeeze(e.raw.t(i, j, k, :));
            L = length(t);

            % Zero mean
            V_shaft = V_shaft - mean(V_shaft);

            Y = fft(V_shaft);
            P2 = abs(Y/L);
            P1 = P2(1:L/2+1);
            P1(2:end-1) = 2*P1(2:end-1);

            f = f_sample * (0:(L/2))/L;
            [~, i_peak] = max(P1);
            f_shaft(i, j, k) = f(i_peak);

%             figure;
%             plot(f, P1)
%             title(['PWM Speed = ' num2str(pwm_speed(k))])
        end
    end
end

% Calculate load cell forces
load_calib = e.calib.load;

V_long_avg = mean(e.raw.V_long, 4);
V_lat_avg = mean(e.raw.V_long, 4);

F_long = load_calib(1, 1) * V_long_avg + load_calib(1, 2);
F_lat = load_calib(2, 1) * V_lat_avg + load_calib(2, 2);

%% Fill new e struct
e_new.geom = geom;
e_new.atm = atm;

e_new.V_supply = V_supply;
e_new.I_supply = I_supply;
e_new.P_total = P_total;
e_new.P_static = P_static;
e_new.T_motor = T_motor;
e_new.f_shaft = f_shaft;
e_new.F_long = F_long;
e_new.F_lat = F_lat;

e_new.aoa = aoa;
e_new.tunnel_switch = tunnel_switch;
e_new.pwm_speed = pwm_speed;

%% Save e_new into new file as e
e = e_new;
save(filename + "_processed.mat", "e");

