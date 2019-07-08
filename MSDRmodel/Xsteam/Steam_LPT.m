pres_inp = [0.00689476, 0.01723689, 0.0344738, 0.05171068, 0.0689476, ...
            0.137895, 0.344738, 1.37895]; % pressure in bar
pres_table = pres_inp*1e5; % onvert to Pascal

m = length(pres_table);

% matrix of zeros to store result
H_f_table = zeros(m); % Enthalpy of sat water
H_v_table = zeros(m); % Enthalpy of sat vapor
H_fg_table= zeros(m); % Enthalpy of vaporization


for i=1:m
    H_f_table(i) = XSteam('hL_p', pres_inp(i))/1e3; % in MJ/kg
    H_v_table(i) = XSteam('hV_p', pres_inp(i))/1e3; % in MJ/kg
    H_fg_table(i) = H_v_table(i) - H_f_table(i);
end