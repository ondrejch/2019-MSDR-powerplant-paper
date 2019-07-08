temp_table = [275, 300, 325, 350, 375, 400, 425, 450, 475, 500, 525, 550, 575, 600, 625, 650, 675, 700]; % temp in deg-C
pres_inp = [30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 200]; % pressure in bar
pres_table = pres_inp*1e-1; % convert to MPa

m = length(temp_table);
n = length(pres_table);

% matrix of zeros to store result
H_s_table = zeros(m, n); % Enthalpy of steam

for i=1:m
    for j=1:n
        H_s_table(i, j) = XSteam('h_pT', pres_inp(j), temp_table(i))/1e3; % in MJ/kg
    end
end

% % Correction for error in the data in XSteam
% H_s_table(1 ,5) = H_s_table(1 ,5) + 1.4;
% H_s_table(1 ,6) = H_s_table(1 ,6) + 1.3;
