% density of steam in the nozzle chest in kg/m^3
rho_table = [10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 100, 125, 150];

% Enthalpy of steam in the nozzle chest
% in kJ/kg
H_inp = [1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400];
H_table = H_inp/1e3; % from kJ/kg to MJ/kg

% length of each array
m = length(rho_table);
n = length(H_inp);

% matrix of zeros to store result
P_nc_table = zeros(n, m);

for i=1:m
    for j=1:n
        % Call XSteam function 'p_hrho' for each value of enthalpy
        % while holding density constant
        % Save the result in each row of the result matrix
        P_nc_table(j,i) = XSteam('p_hrho', H_inp(j), rho_table(i))*1e-1; % bar to MPa
    end
end
    
