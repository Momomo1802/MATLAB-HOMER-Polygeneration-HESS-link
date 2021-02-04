function [simulation_state, matlab_simulation_variables] = MatlabDispatch(simulation_parameters, simulation_state, matlab_simulation_variables)
% load testvariables0.mat

%Shortening Generator variables   
min_load = simulation_parameters.generators(1).minimum_load/100 * simulation_parameters.generators(1).rated_capacity;
max_load = simulation_parameters.generators(1).rated_capacity;

%Shortening Battery & Converter variables
rect_efficiency = simulation_parameters.converters(1).rectifier_efficiency;
inv_efficiency = simulation_parameters.converters(1).inverter_efficiency;
max_SOC = 100; % Maximum percentage of batteries, which is 100%(...duh)
power_remaining = 0;


% Creating Biogas Storage variables and adding generator max load
if simulation_state.current_timestep == 0
    matlab_simulation_variables.BioStor.MAX = max_load * (simulation_parameters.timestep_size_in_seconds / 3600) * 168 ; % Max. Biogas Storage can serve 168 hours / 1 week of maximum generator power. Above this value excess biogas has to be removed.
    matlab_simulation_variables.BioStor.UPLIMIT = 0.95 * matlab_simulation_variables.BioStor.MAX; % This is the upper limit where biogas is consumed to prevent the storage from filling over
    matlab_simulation_variables.BioStor.SOS = 0.5 * matlab_simulation_variables.BioStor.UPLIMIT ; %State of Storage: half full at start of simulation
    matlab_simulation_variables.BioStor.MIN = 0.1 * matlab_simulation_variables.BioStor.UPLIMIT ; % Minimum biogas storage at 10%
    matlab_simulation_variables.gen_rating = max_load; %used for naming the excel file at the end of the simulation
end

% Updating Biogas Storage with production value for corresponding timestep from the excel list
matlab_simulation_variables.BioStor.SOS = matlab_simulation_variables.BioStor.SOS + matlab_simulation_variables.biogas_Production(simulation_state.current_timestep + 1,1); % timestep + 1 to have a real integer

%Increase electric load for biogas compression, cleaning, pumps etc. in relation to biogas production (10% parasitic load of biogas energy content)
simulation_state.ac_bus.load_requested = simulation_state.ac_bus.load_requested ;%+ 0.1 * matlab_simulation_variables.biogas_Production(simulation_state.current_timestep + 1,1);
simulation_state.primary_loads.load_requested = simulation_state.ac_bus.load_requested;

%Setting AC_Bus_Load_Operating_capacity equal to load * (1 + reserve requirement)
simulation_state.ac_bus.operating_capacity_requested = (simulation_state.primary_loads(1).load_requested) * ( 0.95 + (simulation_parameters.operating_reserve.timestep_requirement/100));

%DC_Bus not used to serve loads
simulation_state.dc_bus.operating_capacity_requested = 0;
simulation_state.dc_bus.operating_capacity_served = 0;

%Check for PV, use all power available & reduce load by power from pv
if simulation_parameters.has_pv==true
        simulation_state.pvs(1).power_setpoint = simulation_state.pvs(1).power_available;
        simulation_state.converters(1).inverter_power_output = min (simulation_state.pvs(1).power_setpoint * inv_efficiency/100, simulation_parameters.converters(1).inverter_capacity);
        simulation_state.converters(1).inverter_power_input = simulation_state.converters(1).inverter_power_output / (inv_efficiency/100);
        simulation_state.ac_bus.excess_electricity = simulation_state.ac_bus.excess_electricity + simulation_state.pvs(1).power_setpoint * inv_efficiency/100 - simulation_state.converters(1).inverter_power_output; %PV output above converter capacity is lost
%   disp (simulation_state.ac_bus.load_requested);
        Load_after_PV =  simulation_state.ac_bus.load_requested - simulation_state.converters(1).inverter_power_output;
% disp (Load_after_PV);
        Inv_Cap_Remaining = simulation_parameters.converters(1).inverter_capacity - simulation_state.converters(1).inverter_power_output;
else %No PV, no problem
        Load_after_PV =  simulation_state.ac_bus.load_requested;
        simulation_state.ac_bus.unmet_load = 0; 
        Inv_Cap_Remaining = simulation_parameters.converters(1).inverter_capacity;
        simulation_state.pvs(1).power_setpoint = 0;
end

    
if Load_after_PV <= 0       %PV can serve all demand (dsipatch not needed)
    power_remaining = -1 * Load_after_PV; 
    simulation_state.ac_bus.load_served = simulation_state.ac_bus.load_requested;
    
else    %Run dispatch control  
 if matlab_simulation_variables.BioStor.SOS < matlab_simulation_variables.BioStor.MIN %Bio empty
                simulation_state.generators(1).power_setpoint = 0;
     if  simulation_state.batteries(1).state_of_charge_percent > simulation_parameters.batteries(1).minimum_state_of_charge %Bat empty
         
            if Load_after_PV < min(simulation_state.batteries(1).max_discharge_power *(inv_efficiency/100), Inv_Cap_Remaining)  % Battery can serve full demand
                simulation_state.converters(1).inverter_power_output = simulation_state.converters(1).inverter_power_output + Load_after_PV;
                simulation_state.converters(1).inverter_power_input = (simulation_state.converters(1).inverter_power_output)/(inv_efficiency/100);
                simulation_state.batteries(1).power_setpoint = (-1)*Load_after_PV / (inv_efficiency/100);
                                
            else %Battery cannot serve full demand
                simulation_state.converters(1).inverter_power_output = min(simulation_state.batteries(1).max_discharge_power * inv_efficiency/100, Inv_Cap_Remaining);
                simulation_state.converters(1).inverter_power_input = (simulation_state.converters(1).inverter_power_output)/(inv_efficiency/100);
                simulation_state.batteries(1).power_setpoint = (-1) * min(simulation_state.batteries(1).max_discharge_power, Inv_Cap_Remaining);
            end
     end
     
 elseif matlab_simulation_variables.BioStor.SOS > matlab_simulation_variables.BioStor.MIN && matlab_simulation_variables.BioStor.SOS <= matlab_simulation_variables.BioStor.UPLIMIT % Bio not empty, but also not full

    if simulation_state.batteries(1).state_of_charge_percent < simulation_parameters.batteries(1).minimum_state_of_charge %Battery empty
            simulation_state.generators(1).power_setpoint = max_load ;
             
    elseif simulation_state.batteries(1).state_of_charge_percent >= simulation_parameters.batteries(1).minimum_state_of_charge && simulation_state.batteries(1).state_of_charge_percent < max_SOC %Battery not empty, but also not full     
                            
        if Load_after_PV < min(simulation_state.batteries(1).max_discharge_power * (inv_efficiency/100), Inv_Cap_Remaining)  % Battery can serve full demand
               %disp('battery serves');
                simulation_state.generators(1).power_setpoint = 0;
                simulation_state.converters(1).inverter_power_output = simulation_state.converters(1).inverter_power_output + Load_after_PV;
                simulation_state.converters(1).inverter_power_input = (simulation_state.converters(1).inverter_power_output)/(inv_efficiency/100);
                simulation_state.batteries(1).power_setpoint = (-1)*Load_after_PV / (inv_efficiency/100);
                
         elseif Load_after_PV <= max_load % Generator can serve load and charge battery
                    %disp(Load_after_PV)    ;
                    %disp (min(simulation_state.batteries(1).max_discharge_power * (inv_efficiency/100), Inv_Cap_Remaining));
                    
             %disp('NONONO');
                simulation_state.generators(1).power_setpoint = max_load ;
         elseif Load_after_PV <= max_load + min(simulation_state.batteries(1).max_discharge_power* (inv_efficiency/100), Inv_Cap_Remaining) %Generator with help from Bat can serve load
                 %disp('gen+bat serve')  
             simulation_state.generators(1).power_setpoint = max_load ;
                simulation_state.converters(1).inverter_power_output = simulation_state.converters(1).inverter_power_output + (Load_after_PV - max_load);
                simulation_state.converters(1).inverter_power_input = (simulation_state.converters(1).inverter_power_output)/(inv_efficiency/100);
                simulation_state.batteries(1).power_setpoint = (-1)*(Load_after_PV - max_load)/(inv_efficiency/100);
        else %Demand greater than Generator & Battery together
                simulation_state.generators(1).power_setpoint = max_load ;
                simulation_state.converters(1).inverter_power_output = simulation_state.converters(1).inverter_power_output + min(simulation_state.batteries(1).max_discharge_power * inv_efficiency/100, Inv_Cap_Remaining);
                simulation_state.converters(1).inverter_power_input = (simulation_state.converters(1).inverter_power_output)/(inv_efficiency/100);
                simulation_state.batteries(1).power_setpoint = (-1)* min(simulation_state.batteries(1).max_discharge_power, Inv_Cap_Remaining);
         end
          
    else %Battery full
                         
         if Load_after_PV < min(simulation_state.batteries(1).max_discharge_power * (inv_efficiency/100), Inv_Cap_Remaining)  % Battery can serve full demand
                simulation_state.generators(1).power_setpoint = 0;
                simulation_state.converters(1).inverter_power_output = simulation_state.converters(1).inverter_power_output + Load_after_PV;
                simulation_state.converters(1).inverter_power_input = (simulation_state.converters(1).inverter_power_output)/(inv_efficiency/100);
                simulation_state.batteries(1).power_setpoint = (-1)*Load_after_PV / (inv_efficiency/100);
         elseif Load_after_PV <= max_load % Generator can serve load, but battery cannot
                simulation_state.generators(1).power_setpoint = Load_after_PV ;
         elseif Load_after_PV <= max_load + min(simulation_state.batteries(1).max_discharge_power * (inv_efficiency/100), Inv_Cap_Remaining) %Generator with help from Bat can serve load
                simulation_state.generators(1).power_setpoint = max_load ;
                simulation_state.converters(1).inverter_power_output = simulation_state.converters(1).inverter_power_output + Load_after_PV - max_load;
                simulation_state.converters(1).inverter_power_input = (simulation_state.converters(1).inverter_power_output)/(inv_efficiency/100);
                simulation_state.batteries(1).power_setpoint = (-1)*(Load_after_PV - max_load)/(inv_efficiency/100);
         else %Demand greater than Generator & Battery together
                simulation_state.generators(1).power_setpoint = max_load ;
                simulation_state.converters(1).inverter_power_output = simulation_state.converters(1).inverter_power_output + min(simulation_state.batteries(1).max_discharge_power * inv_efficiency/100, Inv_Cap_Remaining);
                simulation_state.converters(1).inverter_power_input = (simulation_state.converters(1).inverter_power_output)/(inv_efficiency/100);
                simulation_state.batteries(1).power_setpoint = (-1)* min(simulation_state.batteries(1).max_discharge_power, Inv_Cap_Remaining);
         end
    end
%     
 else %Bio full
    
   if simulation_state.batteries(1).state_of_charge_percent < simulation_parameters.batteries(1).minimum_state_of_charge % Battery empty
                    simulation_state.generators(1).power_setpoint = max_load ; 
     
   elseif simulation_state.batteries(1).state_of_charge_percent >= simulation_parameters.batteries(1).minimum_state_of_charge && simulation_state.batteries(1).state_of_charge_percent < max_SOC %Battery not empty, but also not full
                    simulation_state.generators(1).power_setpoint = max_load ;
       
         if         Load_after_PV > max_load % Demand greater than generator can serve alone
                    simulation_state.converters(1).inverter_power_output = simulation_state.converters(1).inverter_power_output + (Load_after_PV - max_load);
                    simulation_state.converters(1).inverter_power_input = (simulation_state.converters(1).inverter_power_output)/(inv_efficiency/100);
                    simulation_state.batteries(1).power_setpoint = (-1)*(Load_after_PV - max_load)/(inv_efficiency/100);
         else       %Demand greater than Generator & Battery together
                    simulation_state.converters(1).inverter_power_output = simulation_state.converters(1).inverter_power_output + min(simulation_state.batteries(1).max_discharge_power * inv_efficiency/100, Inv_Cap_Remaining);
                    simulation_state.converters(1).inverter_power_input = (simulation_state.converters(1).inverter_power_output)/(inv_efficiency/100);
                    simulation_state.batteries(1).power_setpoint = (-1)* min(simulation_state.batteries(1).max_discharge_power, Inv_Cap_Remaining);
         end
        
   else %Battery full
                        
        if         Load_after_PV < min_load % Generator must serve more load than requested
                   simulation_state.generators(1).power_setpoint = min_load ;
        elseif     Load_after_PV <= max_load  % Generator can serve load exactly
                   simulation_state.generators(1).power_setpoint = Load_after_PV ;
        elseif     Load_after_PV <= max_load + min(simulation_state.batteries(1).max_discharge_power * (inv_efficiency/100), Inv_Cap_Remaining) %Generator with help from Bat can serve load
                   simulation_state.generators(1).power_setpoint = max_load ;
                   simulation_state.converters(1).inverter_power_output = simulation_state.converters(1).inverter_power_output + (Load_after_PV - max_load);
                   simulation_state.converters(1).inverter_power_input = (simulation_state.converters(1).inverter_power_output)/(inv_efficiency/100);
                   simulation_state.batteries(1).power_setpoint = (-1)*(Load_after_PV - max_load)/(inv_efficiency/100);
        else %Demand greater than Generator & Battery together
                   simulation_state.generators(1).power_setpoint = max_load ;
                   simulation_state.converters(1).inverter_power_output = simulation_state.converters(1).inverter_power_output + min(simulation_state.batteries(1).max_discharge_power * inv_efficiency/100, Inv_Cap_Remaining);
                   simulation_state.converters(1).inverter_power_input = (simulation_state.converters(1).inverter_power_output)/(inv_efficiency/100);
                   simulation_state.batteries(1).power_setpoint = (-1)* min(simulation_state.batteries(1).max_discharge_power, Inv_Cap_Remaining);
        end
   end
 end

 
 %Setting total ac_bus load served to Bat + Eng + PV
 simulation_state.ac_bus.load_served = simulation_state.converters(1).inverter_power_output + simulation_state.generators(1).power_setpoint;
 

 
 %Set unmet load if not enough electricity generated
 if simulation_state.ac_bus.load_requested > simulation_state.ac_bus.load_served
        simulation_state.ac_bus.unmet_load = simulation_state.ac_bus.load_requested - simulation_state.ac_bus.load_served;
 else
        simulation_state.ac_bus.unmet_load = 0;  
        power_remaining = simulation_state.ac_bus.load_served - simulation_state.ac_bus.load_requested; 
        simulation_state.ac_bus.load_served = simulation_state.ac_bus.load_requested;
 end
end 
 
 %Equaling Primary loads to AC_load
 simulation_state.primary_loads(1).load_served = simulation_state.ac_bus.load_served; 
  
%Use remaining power to charge battery or, if unavailabe, add to excess 
if power_remaining > 0 && simulation_state.batteries(1).state_of_charge_percent < max_SOC
      if power_remaining > min(simulation_state.batteries(1).max_charge_power /(rect_efficiency/100), simulation_parameters.converters(1).rectifier_capacity) %More power remaining than can be used for charging battery
        simulation_state.converters(1).rectifier_power_input = min((simulation_state.batteries(1).max_charge_power /(rect_efficiency/100)), simulation_parameters.converters(1).rectifier_capacity);
        simulation_state.converters(1).rectifier_power_output = simulation_state.converters(1).rectifier_power_input * (rect_efficiency/100);
        simulation_state.batteries(1).power_setpoint = simulation_state.batteries(1).power_setpoint + simulation_state.converters(1).rectifier_power_output;
        simulation_state.ac_bus.excess_electricity = simulation_state.ac_bus.excess_electricity + power_remaining - simulation_state.converters(1).rectifier_power_input; 
        
      else %all remaining power can be used for charging battery
        simulation_state.converters(1).rectifier_power_input = power_remaining;
        simulation_state.converters(1).rectifier_power_output = simulation_state.converters(1).rectifier_power_input * (rect_efficiency/100);
        simulation_state.batteries(1).power_setpoint = simulation_state.batteries(1).power_setpoint + simulation_state.converters(1).rectifier_power_output;
      end
else %Battery full or no remaining power
    simulation_state.ac_bus.excess_electricity = simulation_state.ac_bus.excess_electricity + power_remaining;
end


% Setting operating capacity using pv if available
if simulation_parameters.has_pv==true
    simulation_state.ac_bus.operating_capacity_served = min (simulation_state.pvs(1).power_setpoint * inv_efficiency/100, simulation_parameters.converters(1).inverter_capacity);
end

% Adding operating capacity using engine if available
if matlab_simulation_variables.BioStor.SOS > matlab_simulation_variables.BioStor.MIN
    simulation_state.ac_bus.operating_capacity_served = simulation_state.ac_bus.operating_capacity_served + max_load;
end

% Adding operating capacity using battery if available
if simulation_state.batteries(1).state_of_charge_percent >= simulation_parameters.batteries(1).minimum_state_of_charge
    simulation_state.ac_bus.operating_capacity_served = simulation_state.ac_bus.operating_capacity_served + min(simulation_state.batteries(1).max_discharge_power * inv_efficiency/100, Inv_Cap_Remaining);
end

% Setting if enough capacity is available or not
if simulation_state.ac_bus.operating_capacity_served >= simulation_state.ac_bus.operating_capacity_requested
     simulation_state.ac_bus.operating_capacity_served = simulation_state.ac_bus.operating_capacity_requested;
     simulation_state.ac_bus.capacity_shortage = 0;
else 
     simulation_state.ac_bus.capacity_shortage = simulation_state.ac_bus.operating_capacity_requested - simulation_state.ac_bus.operating_capacity_served;
end
 
% Reducing Biogas storage  by fuel consumption using HOMER Slope and Intercept Coeff.
Interc_Co = simulation_parameters.generators(1).fuel_curve_intercept; %[kg/hr/kW rated]
Slope = simulation_parameters.generators(1).fuel_curve_slope; %[kg/hr/kW output]

if simulation_state.generators(1).power_setpoint ~= 0
Fuel_usage = (Interc_Co * max_load + Slope * simulation_state.generators(1).power_setpoint) * matlab_simulation_variables.biogas_energy_density ;%LHV (kWh/kg) of Biogas from Excel sheet
else 
Fuel_usage = 0;
end
matlab_simulation_variables.BioStor.SOS = matlab_simulation_variables.BioStor.SOS - Fuel_usage ;

% Calculating useful thermal exhaust heat
if simulation_state.generators(1).power_setpoint ~= 0
Thermal_energy = (Fuel_usage - (simulation_state.generators(1).power_setpoint * simulation_parameters.timestep_size_in_seconds/3600))*0.7; %Assuming heat recovery is 70%
else
Thermal_energy = 0;
end

%For the water sales to appear in the economics electric sales are used in HOMER instead
PW_Price_per_m3 = 145; %Potable water price per volume in $/m3
simulation_state.grids(1).sellback_rate = PW_Price_per_m3; %workaround for water sales
E_need_WaterPuri = 800; %Thermal energy need of the water purifier in kWh/m3
matlab_simulation_variables.WaterProd = 0;
matlab_simulation_variables.WaterProd = Thermal_energy/E_need_WaterPuri; %WaterProduction in m3
simulation_state.grids(1).grid_sales = matlab_simulation_variables.WaterProd; %workaround for water production

%Calculating biogas excess if there is and selling everything to Upper Limit
if matlab_simulation_variables.BioStor.SOS > matlab_simulation_variables.BioStor.MAX
matlab_simulation_variables.BioStor.Sales = matlab_simulation_variables.BioStor.UPLIMIT - matlab_simulation_variables.BioStor.SOS;
matlab_simulation_variables.BioStor.SOS = matlab_simulation_variables.BioStor.UPLIMIT;
else
matlab_simulation_variables.BioStor.Sales = 0;
end

 
end



