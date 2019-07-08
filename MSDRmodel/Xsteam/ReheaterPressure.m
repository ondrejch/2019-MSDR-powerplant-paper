% density of steam in the reheater in kg/m^3
rho_rh_table = [1, 2.5, 5, 7.5, 10, 15, 20];

% Enthalpy of steam in the reheater
H_rh_inp = [100, 200, 500, 750, 1000, 1500, 2000, 2500, 3000, 3500];
H_rh_table2 = H_rh_inp./1e3; % from kJ/kg to MJ/kg

% length of each array
m = length(rho_rh_table);
n = length(H_rh_inp);

% matrix of zeros to store result
P_rh_table = zeros(n, m);

for i=1:m
    for j=1:n
        % Call XSteam function 'p_hrho' for each value of enthalpy
        % while holding density constant
        % Save the result in each row of the result matrix
        P_rh_table(j,i) = XSteam('p_hrho', H_rh_inp(j), rho_rh_table(i))*1e-1; % bar to MPa
    end
end
    
