% Pressure of saturated steam in MPa
P_rh_table1 = [0.01, 0.02, 0.05, 0.1, 0.25, 0.5, 1, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 10.0, 12.5, 15.0];

% length of each array
m = length(P_rh_table1);

% matrix of zeros to store result
SV_nc_table = zeros(1,m);

for i=1:m
    % Call XSteam function 'vV_p' for each value of enthalpy
    % while holding density constant
    % Save the result in each row of the result matrix
    SV_nc_table(1,i) = XSteam('vV_p', P_rh_table1(i).*1e1); % Pressure from MPa to bar
end
    
