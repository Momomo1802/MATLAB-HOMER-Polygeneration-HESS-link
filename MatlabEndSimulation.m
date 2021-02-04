function myErrs = MatlabEndSimulation(simulation_parameters, matlab_simulation_variables)
myErrs.simulation_errors = {};
myErrs.simulation_warnings = {};

% Adds more lines in the result.xls for yearly and daily sums
matlab_simulation_variables.mat = [matlab_simulation_variables.mat; "Total yearly sums:", "=sum(B2:B8761)", "=sum(C2:C8761)", "=sum(D2:D8761)", "=sum(E2:E8761)", "=sum(F2:F8761)", "=sum(G2:G8761)"];
matlab_simulation_variables.mat = [matlab_simulation_variables.mat;"Daily values", "=B8762/365", "=C8762/365",  "=D8762/365",  "=E8762/365",  "=F8762/365",  "=G8762/365"];  

%Cost capital has been used for battery value in excel naming
Bat.trick = round(simulation_parameters.batteries(1).cost.capital);

% Optional: Writes an excel file named after the component sizes
% xlswrite(['Gen',num2str(matlab_simulation_variables.gen_rating), '&PV',num2str(simulation_parameters.pvs(1).rated_capacity), '&Bat', num2str(Bat.trick), '&Conv', num2str(simulation_parameters.converters(1).inverter_capacity)], matlab_simulation_variables.mat);

end
