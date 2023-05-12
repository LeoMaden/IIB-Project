
filename = "run2";

data = load(filename + ".mat");
e = data.e;
N = data.N;

out = e;

out = rmfield(out, "long_force");
out = rmfield(out, "lat_force");
out = rmfield(out, "J");
out = rmfield(out, "calib");
out = rmfield(out, "rho");
out = rmfield(out, "env");

% FFT of shaft freq to get average


out.raw = rmfield(out.raw, "V_shaft");
out.raw = rmfield(out.raw, "t");
out.raw = rmfield(out.raw, "V_long");
out.raw = rmfield(out.raw, "V_lat");
out.raw = rmfield(out.raw, "V");
out.raw = rmfield(out.raw, "pxie_time");
out.raw = rmfield(out.raw, "V_Long_av");
out.raw = rmfield(out.raw, "V_Lat_av");

out.raw.V_long_avg = e.raw.V_Long_av;
out.raw.V_lat_avg = e.raw.V_Lat_av;

out.load_calib = e.calib.load;


% Fix run1
if filename == "run1"
    out.aoa(end-1:end) = [];
    out.tunnel_switch(end) = [];

    out.raw.V_long_avg(2:3, :, :) = [];
    out.raw.V_long_avg(:, end, :) = [];

    out.raw.V_lat_avg(2:3, :, :) = [];
    out.raw.V_lat_avg(:, end, :) = [];

    out.raw.P(2:3, :, :, :) = [];
    out.raw.P(:, end, :, :) = [];
    out.raw.P(:, :, :, 3:end) = [];

    out.raw.T_motor(2:3, :, :) = [];
    out.raw.T_motor(:, end, :) = [];

    out.raw.I_supply(2:3, :, :) = [];
    out.raw.I_supply(:, end, :) = [];

    out.raw.shaft_freq(2:3, :, :) = [];
    out.raw.shaft_freq(:, end, :) = [];
end

e = out;

save(filename + "_proc.mat", "e", "N");

