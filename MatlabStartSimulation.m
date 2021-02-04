function [myErr, matlab_simulation_variables] = MatlabStartSimulation(simulation_parameters)
myErr.error_description = '';
myErr.severity_code = '';

% Reads in hourly biogas production values and the energy density from Excel
sheet = 1;
xlRange1 = 'A2:A8761';
xlRange2 = 'E11:E11';
matlab_simulation_variables.biogas_Production = xlsread('BiogasProductionCalculator.xls',sheet,xlRange1);
matlab_simulation_variables.biogas_energy_density = xlsread('BiogasProductionCalculator.xls',sheet,xlRange2);

%Creating Biostorage variables and exhaust heat variable
matlab_simulation_variables.BioStor.Max = [];
matlab_simulation_variables.BioStor.UPLIMIT = [];
matlab_simulation_variables.BioStor.SOS = [];
matlab_simulation_variables.BioStor.MIN = [];
matlab_simulation_variables.Exhaustheat = [] ;

%Creating Water variables
matlab_simulation_variables.WaterProd = [];

%Creating variables for excel file and naming columns
matlab_simulation_variables.mat = ["Timestep","Biostor SOS (kWh)", "Fuel usage (kWh)", "Biogas Sales (kWh)", "Power remaining (kWh)", "Thermal Energy (kWh)", "Water production(m3)"];

end