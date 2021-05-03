%------------------- Panagiotis Manias ------------------------------------
clc
clear all
close all
    %-------------------------Default Input Data---------------------
%----------Fuel Cell inputs------
Fuel_utilization=0.9  %the SOFC efficiency
T_operating =1000      %Operating temperature (C)
n_con=0.95            %AC/DC Converter Efficiency
P_Spec=54             %Specified Power Requirement (MW)
b=0.025               %Electrodes activation potential constant
m=0.00001             %Concentration loss constant
n=0.15                %Concentration loss index
n_el=2                %Number of electrons exchange in the full REDOX process
Res=0.0005            %Ohmic Resistance (Ω)
i=15                  %Cell Current (A)
V_Ac=b*log(i)         %Activation Voltage    (V)
V_con=m*exp(n*i)      %Concentration Voltage (V)
I_s=4000              %Amperage supply per stack (A)
n_WHRS=0.25            %Heat Recovery Effectiveness

%---------- Solid Oxide Fuel Cell Parameters ---------------------

T_ref=273+T_operating         %Operating temperature (K)
DG=-170                       % Average Gibbs Free energy change (KJ/mol)at 800 C %192 % 181 % 170
R=8.314/1000;                 %Particular Gas Constant of Air                (kJ/mol.K)
K_p=1;                        %Product of partial pressures
K_o=exp(DG/(R*T_ref));        %Equilibrium Constant
DH=-251;                      %Enthalpy Change  of water shift reaction at reference temp 800 C (kJ/mole)%249 %250 % 251
F=96497;                      %Fahrad Number (coulombs/g-mole electron)
FCS=500;                      %Number of fuel cells per stack
t=0.12                        %Fuel Cell Thickness (m)
A_c=0.196                     %Fuel Cell Area      (m^2)
density=3000;                 %Solid Oxide Fuel Cell Density (kg/m^3)

%---------Fuel Inputs/Combustion Efficiency----

MR_H2O=18                                  %gr/mol
m_H2=I_s*1000/(n_el*F*MR_H2O)                %Theoretical Hydrogen fuel input kg/s
MR_CH4=16.4                   %gr/mol
MR_CO2=44.0                    %gr/mol

%-----------------------------------------------------------------

W_out=(R*T_ref*log(K_o)-R*T_ref*log(K_p)); %Specific Mole Work output (kJ/mole)

W=W_out*1000/MR_H2O;                       %Brake specific Work per stack    (kJ/kg)

%------------------------------ Overal Fuel Cell Performance----------
Theoretical_eff= 100*W_out/DH;                            %Theoretical Efficiency
W_b=-W;                                                   %Net Brake specific Work     (kJ/kg)
Voc=-W_out*1000/(n_el*F);                                 %Open Circuit Voltage (V)
V_Cell=Voc-Res*i-V_Ac-V_con;                              %Cell Voltage       (V)
n_v=V_Cell/Voc;                                           %Voltage Efficiency                                          %Current Supply
P_stack=V_Cell*I_s;                                       %Power From each stack of Fuel cell (Kw)
P_S_N=(P_Spec*1000)/(P_stack);                              %Number of Stacks
Cell_Volume=P_S_N*t*A_c*FCS;                              %Total Volume of Cells (m^3)
P_Ttot=P_S_N*P_stack;                                     %Total Power Output (KW)
W_hdat=(W_b-(DH*1000/MR_H2O))*m_H2/1000*Fuel_utilization; %Heat Loss          (MW)
M_fdat=m_H2/Fuel_utilization;                             %Rate of Fuel Consumption(kg/sec)
bsfc=(M_fdat*1000/(P_Ttot))*3600;                         %Brake Specific Fuel Consumption                      (gr/KW-hr)
MFD=M_fdat*3600*24/1000;                                  %Mass of Fuel consumed per day  (Tonnes)
V_stack=P_stack/I_s;                                      %Voltage supply per Stack (V)
Actual_efficiency= (1-(W_hdat*1000/P_Ttot))*100*Fuel_utilization;          %Actual efficiency value of the fuel cell system (%)
Weight=Cell_Volume*density/1000;                          %Total system Weight (Tonnes)
CO2_Emission=bsfc*MR_CO2/MR_CH4;                          %CO2 emissions 

%--------------------With WHRS--------------------------

P_WHRS=n_WHRS*W_hdat*1000;                                      %Wasted heat recovered (KW)
W_b2=W_b+W_b*n_WHRS;
P_Ttot2=P_Ttot-P_WHRS;                                          %Total Power Output (KW)
M_fdat2=((P_Ttot2)/P_Ttot)*m_H2/Fuel_utilization;                       %Rate of Fuel Consumption(kg/sec)
bsfc2=(M_fdat2*1000/(P_Ttot2))*3600;                            %Brake Specific Fuel Consumption                      (gr/KW-hr)
MFD2=M_fdat2*3600*24/1000;                                      %Mass of Fuel consumed per day (Tonnes)
n_WHRS=(1-((W_hdat*1000-P_WHRS)/(P_Ttot)))*100*Fuel_utilization;
%--------------------Display Outputs--------------------------
disp('                                                                   ');
disp('                                                                   ');

disp('---------------------------Fuel Properties-------------------------- ');


%HV=['The Fuel Energy value is: ',num2str(FHV), ' kJ/kg'];
%disp(HV)


disp('                                                                   ');
disp('                                                                   ');

disp('---------------------------Cycle Operating Point Properties-------------------------- ');
disp('                                                                   ');
disp('                                                                   ');
aa=['The Maximum Temprature of The Cycle is: ',num2str(T_ref), ' Deg K'];
    disp(aa)


disp('                                                                   ');
disp('                                                                   ');

disp('---------------------------Overal Fuel Cell Performances-------------------------- ');
disp('                                                                   ');
disp('                                                                   ');

a0=['The output is: ',num2str(W_out), ' kJ/mol'];
    disp(a0)

a1=['The Theoretical Efficiency is: ',num2str(Theoretical_eff), ' %'];
    disp(a1)
    
a15=['The Actual Efficiency is: ',num2str(Actual_efficiency), ' %'];
    disp(a15)

a2 = ['The open circuit voltage is: ',num2str(Voc), ' V'];
    disp(a2)

a3 = ['The Cell voltage is: ',num2str(V_Cell), ' V'];
    disp(a3)
    
a7 = ['The Power provided by the system is: ',num2str(P_Ttot), ' KW'];
    disp(a7)
    
a6 = ['The power supply per stack is: ',num2str(P_stack), ' kW'];
    disp(a6)   
    
a9 = ['The heat generated from the system is: ',num2str(W_hdat), ' MW'];
    disp(a9)

a10 = ['The total volume of the fuel stacks combined is: ',num2str(Cell_Volume), ' m^3'];
    disp(a10)   
    
a11 = ['The total number of stacks is: ',num2str(P_S_N)];
    disp(a11)   
    
a4=['The Fuel Consumption rate is : ',num2str(M_fdat), ' kg/sec'];
    disp(a4)
    
a8=['The Brake Specific Fuel Consumption bsfc is : ',num2str(bsfc), ' gr/KW-hr'];
    disp(a8)

a5=['The work output per kg of H2 fuel is : ',num2str(W_b), 'KJ/kg'];
    disp(a5)

a12=['The Fuel Consumption per day is : ',num2str(MFD), ' Tonnes'];
    disp(a12)    
    
a13=['The Voltage Per Stack is : ',num2str(V_stack), ' V'];
    disp(a13)  
    
a14=['The total Amperage is : ',num2str(I_s), ' A'];
    disp(a14) 
    
a20=['The total weight of the system is : ',num2str(Weight), ' Tonnes'];
    disp(a20) 

a21 = ['The CO2 emissions per kWh are : ',num2str(CO2_Emission), ' gr.'];
    disp(a21)
    
disp('---------------------------Overal Fuel Cell Performance with WHRS-------------------------- ');
disp('                                                                   ');
disp('                                                                   ');

a16 = ['The heat recovered from the system is: ',num2str(P_WHRS), ' KW'];
    disp(a16)

a17=['The new Fuel Consumption rate is : ',num2str(M_fdat2), ' kg/sec'];
    disp(a17)

a18=['The new Fuel Consumption per day is : ',num2str(MFD2), ' Tonnes'];
    disp(a18)    

a19=['The new system efficiency is : ',num2str(n_WHRS), ' %'];
    disp(a19)   
%-------------------------Display Statements-------------------------

