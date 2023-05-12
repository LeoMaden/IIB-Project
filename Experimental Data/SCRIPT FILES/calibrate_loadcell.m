%% Calibrate load Cells

% The script measures the voltage for zero load, and for a test load case,
% and produces a matrix of the form
%
% M_longitudinal    C_longitudinal
% M_lateral         C_lateral


loadV = [0 0; 0 0];         % matrix for experimental voltage readings
load_calib = [0 0; 0 0];    % output matrix
calib_mass = 1;             % Test Mass in kg

% Read the zero load voltages in the load cells
uiwait(msgbox('Remove all calibration loads, align fan with freestream direction'))

[q,t] = exp_pxie_read(s.pxie);
loadV(1,1) = mean(q(:,N.V_long));
loadV(2,1) = mean(q(:,N.V_lat));

% Read the voltage with -9.81*calib_mass on the longitudinal load cell
uiwait(msgbox(['Add calibration mass of ' num2str(calib_mass) ' kg']));
[q,t] = exp_pxie_read(s.pxie);
loadV(1,2) = mean(q(:,N.V_long));

% Read the voltage with -9.81*calib_mass on the lateral load cell
uiwait(msgbox(['Remove Mass, align fan to transverse direction, add calibration mass of ' num2str(calib_mass) ' kg']));
[q,t] = exp_pxie_read(s.pxie);
loadV(2,2) = mean(q(:,N.V_lat));

load_calib(1,1) = (0 + 9.81*calib_mass)./(loadV(1,1) - loadV(1,2));
load_calib(1,2) = - load_calib(1,1) * loadV(1,1);

load_calib(2,1) = (0 + 9.81*calib_mass)./(loadV(2,1) - loadV(2,2));
load_calib(2,2) = - load_calib(2,1) * loadV(2,1);


save('load_calib.mat', 'load_calib','calib_mass');