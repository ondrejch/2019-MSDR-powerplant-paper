% Pressure of saturated region in the OTSG in MPa
P_sat_table = [0.01, 0.25, 0.5, 1, 2.5, 5, 7.5, 10, 12.5, 15, 20, 22];
P_sat_inp = [0.01, 0.25, 0.5, 1, 2.5, 5, 7.5, 10, 12.5, 15, 20, 22].*1e1;

% length of each array
m = length(P_sat_table);

% matrix of zeros to store result
T_sat_table = zeros(m,1);

for i=1:m
        % Call XSteam function 'Tsat_p' for each value of Sat Pressure
        % Save the result in each row of the result matrix
        T_sat_table(i,1) = XSteam('Tsat_p', P_sat_inp(i));
end
    
