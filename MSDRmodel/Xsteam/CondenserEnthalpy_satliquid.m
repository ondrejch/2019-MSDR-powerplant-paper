% Pressure of saturated steam in MPa
P_con_table = [1e-3, 2.5e-3, 5e-3, 7.5e-3, 1e-2, 5e-2, 1e-1];

% length of each array
m = length(P_con_table);

% matrix of zeros to store result
H_con_table = zeros(1,m);

for i=1:m
    % Call XSteam function 'hL_p' for each value of enthalpy
    % while holding density constant
    % Save the result in each row of the result matrix
    H_con_table(1,i) = XSteam('hL_p', P_con_table(i).*1e1)./1e3; % Pressure from MPa to bar, Enthalpy from kJ/kg to MJ/kg
end
    
