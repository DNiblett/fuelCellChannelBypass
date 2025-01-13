%% 0D Fuel Cell Model %%
% Collaboration between Newcastle University and University of New South Wales (Dr Quentin Meyer)
% Developed by Dr Daniel Niblett 2023

% Assumptions: Steady-State, 0D, Half rib, Half Channel Unit Cell, Constant Properties
% 1) Water Condensation as film thickness under flow field rib - Changes effective diffusivity
% 2) Electrolyte-Electrode Charge Transfer - Butler-Volmer
% 3) Ohm's Law for Electrolyte (membrane) conduction
% 4) Mass transport accounted for by adjusting exchange current density from ambient conditions
% 5) Constant Current is applied
% 6) Cell voltage is difference between anode and cathode sides

clc 
clear 

plotExperimentalData = 0; % make 1 if you want to plot the extracted data from one of the experiments
tau_list = [1.8;1.8;1.8];
o2conc = [0.21;0.21];
waterThickness = [209e-6;0e-6];
lineCol = {'k';'r';'b'};

% Repeating the polarisation curve twice:
% 1) Water film is removed by by-pass channels
% 2) Water film remains and accumulates under the rib

for repeat = 1:2

alphac = 0.5;                   % Cathode Charge Transfer Coefficient
alphaa = 0.5;                   % Anode Charge Transfer Coefficient
constants.ic = 0.3;             % Exchange Current Density Cathode [A/m2]
constants.ia = 0.3;             % Exchange Current Density Anode [A/m2]
constants.tafelORR = 125./2303; % Tafel Slope [V/dec]
constants.bc = 0.025;           % Cathode Tafel Slope [V]
constants.ba = 0.025;           % Anode Tafel Slope [V]
Asa = 250;                      % Anode Catalyst Roughness [m2/m2]
Asc = 250;                      % Cathode Catalyst Roughness [m2/m2]
constants.temperature = 80+273; % Fuel Cell Operating Temperature [K]
constants.gasConstant = 8.314;  % Ideal Gas Constant [J/[K.mol]]
constants.faraday = 96485;      % Faraday's Constant [C/mol]
constants.Do2 = 2e-5;           % Oxygen Diffusion Coefficient in air [m2/s]
constants.sigma.carbon = 1000;  % GDL Electroic Conductivity [S/m];
constants.sigma.membrane = 13;  % Membrane Ionic Conductivity [S/m];
constants.Erev = 1.23;          % Equalibrium potential [V]
E0a = 0;
E0c = 1.23;

% Geometry 
geo.anode.GDL.length = 199e-6;
geo.anode.GDL.porosity = 0.78;
geo.cathode.GDL.length = 199e-6;
geo.cathode.GDL.porosity = 0.78;
geo.membrane.length = 25.4e-6;
geo.channel.height = 1000e-6;
geo.GDL.thickness = 199e-6;
geo.MPL.thickness = 5e-6;
geo.CL.thickness = 5e-6;
geo.MPL.porosity = 0.3;
geo.CL.porosity = 0.3;
geo.GDL.porosity = 0.78;
tau = 2;

% Calculate Effective Diffusivity Using Half-Channel, Half-Rib Analytical Equation
delta_w = waterThickness(repeat);
delta_g = geo.GDL.thickness;
Wc = geo.channel.height./2;
Wr = geo.channel.height./2;
Dratio = delta_g./((Wc+Wr).* ((delta_w./Wc) + (sqrt((Wc/2).^2 + (delta_g -delta_w).^2)./sqrt((Wc+Wr).^2 ...
     + (delta_g -delta_w).^2))));
D1 = (constants.Do2 + (Dratio.*constants.Do2.*geo.GDL.porosity.^tau))./2;
D2 = ((Dratio.*constants.Do2.*geo.GDL.porosity.^tau) + ((Dratio.*constants.Do2.*geo.MPL.porosity.^tau)))./2;
D3 = ((Dratio.*constants.Do2.*geo.MPL.porosity.^tau) + ((Dratio.*constants.Do2.*geo.CL.porosity.^tau)))./2;
L1 = ((geo.channel.height)/2)+(geo.GDL.thickness./2);
L2 = ((geo.GDL.thickness)/2)+(geo.MPL.thickness./2);
L3 = ((geo.MPL.thickness)/2)+(geo.CL.thickness./2);

% Boundary/Initial Conditions
conditions.pressure = 101325;
conditions.o2fraction = o2conc(repeat);
tau = tau_list(repeat);
%pv=nRT n/v = p/RT
C_o2 = (conditions.pressure.*conditions.o2fraction)./(constants.gasConstant.*constants.temperature);
C_ref = 7.2; % Reference Concentration of Oxygen in air [mol/m3]
t = 0;
count = 0;
for i = 10:10:40000
t=t+1;
eta_guess = 1;
ba = constants.tafelORR;
i0a = constants.ia;
bc = constants.tafelORR;
i0c = constants.ic;
ic = i;

% Calculate oxygen concentration at catalyst layer (cl)
C_cl = C_o2 - i.*(L1./D1 + L2./D2 + L3./D3)./(4.*constants.faraday);
C_ratio  = (C_cl/C_ref);
C_ratio(C_ratio<0)=1e-10;

phi_sc_end = 0;
phi_sc = phi_sc_end - (i*geo.anode.GDL.length/constants.sigma.carbon);
phi_ec = -bc.*log((ic./Asc)/(i0c.*C_ratio)) + phi_sc + E0c;

% Butler-Volmer Cathode
eta_next = eta_guess;
j = ic;
alpha = constants.gasConstant.*constants.temperature ./ (constants.bc*constants.faraday);
alpha = 0.5;
j0 = constants.ic.*Asc.*C_ratio; 
F = constants.faraday;               
R = constants.gasConstant;   
T = constants.temperature;

% Local Charge Balance j_k == j_local -> 0 = j_k - j_local
j_butler_volmer = @(eta) j0 * (exp(alpha * F * eta / (R * T)) - exp(-(1 - alpha) * F * eta / (R * T))) - j;
j_butler_volmer_prime = @(eta) j0 * (alpha * F / (R * T) * exp(alpha * F * eta / (R * T)) + (1 - alpha) * F / (R * T) * exp(-(1 - alpha) * F * eta / (R * T)));

tolerance = 1e-6;
% Maximum number of iterations
max_iterations = 10000;
relaxationFactor = 0.5;
% Newton Method
for jj = 1:max_iterations
    eta_next = eta_guess - (j_butler_volmer(eta_guess) / j_butler_volmer_prime(eta_guess)).*relaxationFactor;
    if abs(eta_next - eta_guess) < tolerance
        break;
    end
    eta_guess = eta_next;
end
phi_ec = -eta_guess + phi_sc + E0c;
phi_ea = phi_ec - (i*geo.membrane.length/constants.sigma.membrane);

% Butler-Volmer Anode

eta_next = eta_guess;
j = ic;
alpha = 0.5;
j0 = constants.ia.*Asa; 
F = constants.faraday;  
R = constants.gasConstant;  
T = constants.temperature;
j_butler_volmer = @(eta) j0 * (exp(alpha * F * eta / (R * T)) - exp(-(1 - alpha) * F * eta / (R * T))) - j;
j_butler_volmer_prime = @(eta) j0 * (alpha * F / (R * T) * exp(alpha * F * eta / (R * T)) + (1 - alpha) * F / (R * T) * exp(-(1 - alpha) * F * eta / (R * T)));
tolerance = 1e-6;
max_iterations = 10000;
relaxationFactor = 0.5;
for jj = 1:max_iterations
    eta_next = eta_guess - (j_butler_volmer(eta_guess) / j_butler_volmer_prime(eta_guess)).*relaxationFactor;
    if abs(eta_next - eta_guess) < tolerance
        break;
    end
    eta_guess = eta_next;
end
phi_sa = phi_ea - E0a - eta_guess;
phi_sa_end = phi_sa - (i*geo.anode.GDL.length/constants.sigma.carbon);
Ucell = phi_sa_end - phi_sc_end;

if Ucell<=0
    count = count + 1;
    if count == 2
        break
    end
end

Ucell(Ucell<0)=0;
sr(t,1) = i;
sr(t,2) = Ucell;
sr(t,3) = C_cl;

end


figure(1)
plot(sr(:,1)./10000,sr(:,2),'linewidth',3,'Color',string(lineCol(repeat)))
hold on
ylabel('Cell Voltage (V)')
xlabel('i (A cm^{-2})')
set(gca,'fontsize',16)


if repeat == 1
C_flood = sr(:,3);
i_flood = sr(:,1);
V_flood = sr(:,2);
end

if repeat == 2
C_dry = sr(:,3);
i_dry = sr(:,1);
V_dry = sr(:,2);
end



end


%% Experimental Data %%

% If plot experimental data is on

if plotExperimentalData == 1
yyaxis left

regularDataAir = [0.009122904770624274, 0.9813767374353639
0.017551565587517914, 0.8510341201805933
0.02973116582432478, 0.8302191586617044
0.05115499501784471, 0.8115914026086947
0.07875287723677438, 0.7962472631847374
0.11556128724961434, 0.7776138902712606
0.1984290764645682, 0.7436296377012661
0.30288246045449385, 0.717304536064176
0.3981253166505089, 0.6942686679166054
0.4964585694755087, 0.6734222519791008
0.6009456546282369, 0.6525735892974094
0.6992856476857973, 0.6328224611509846
0.7976323809759182, 0.6141666207956396
0.8990223292671008, 0.5900332181128023
1.0035094144198289, 0.5691845554311108
1.1018224465471476, 0.545052276120367
1.2001422189070263, 0.522015284600703
1.2984754717320266, 0.5011688686631985
1.3937048474629203, 0.4759424249334684
1.4951015359866635, 0.45290431004171083
1.6118625846320453, 0.426574714916247
1.7040083039665146, 0.4002541067675306
1.7961270623707422, 0.36955234745449506
1.9097572729917758, 0.3344615733724866
2.0018558106983217, 0.30047395068621185
2.0969571220105663, 0.2544370389259666
2.1981920449528576, 0.20511201704829496
2.305459476036788, 0.136069568187001];

minichannelData_Air = [0.012253742794972267, 0.9901379163919085
0.020716104774668342, 0.8652717380925367
0.03285526361611235, 0.8378850498271692
0.06959627130334745, 0.8082987990028949
0.10641142154874816, 0.790760713880498
0.20468401228070376, 0.7600567078232756
0.30298356394290116, 0.7337338529303723
0.40749086979331073, 0.7161710536219201
0.5027539466870071, 0.6964210488475888
0.6011074202096885, 0.6788604962833235
0.7056282065252192, 0.6634882725570308
0.8009182443491573, 0.6481194189470184
0.9023621145008239, 0.6327483185928192
1.0007290684886265, 0.6173783416107134
1.1021796788728535, 0.603102529047594
1.2005601133257766, 0.5899231276476476
1.298927067313579, 0.5745531506655419
1.4034613340942306, 0.5613715025214087
1.5049254249435786, 0.5492862655404487
1.6002222030000777, 0.535012699721516
1.7016728133843046, 0.5207368871583966
1.8000532478372278, 0.5075574857584503
1.8984336822901513, 0.49437808435850406
2.0029544686056817, 0.4790058606322114
2.0982377661970597, 0.4625417192311193
2.199681636348726, 0.4471706188769201
2.3041889421991355, 0.42960781956846794
2.3994789800230745, 0.41423896595845544
2.50091610994218, 0.39777257781317654
2.596172446603316, 0.37692728524776553
2.7006662719886045, 0.35717391035715373
2.799006265046165, 0.33742278221072897
2.906543305432515, 0.3121918449926251
2.9986957649995447, 0.2869665246349884
3.0969885764291822, 0.2595483819510054
3.2014217397214257, 0.22993741694067593
3.3027442856870044, 0.19485113634704132];

minichannelData_o2 = [0.009190307096229144, 0.9923296153461615
0.01502734849361409, 0.9408488424212262
0.014831881749359899, 0.9090854964799135
0.014616194307424224, 0.8740362871653615
0.029778347452248277, 0.8378861731992626
0.06655305630228592, 0.8137763613303871
0.11869886550652275, 0.787470357018885
0.1616678480796515, 0.7699300251523011
0.2353925118262995, 0.7501878839826237
0.39521690630065687, 0.7216519860656925
0.508988661805461, 0.709562255596359
0.6012220041632169, 0.6974803887316792
0.7057562709438684, 0.684298740587546
0.8072271020257767, 0.6733087913976659
0.9086979331076853, 0.6623188422077857
1.013238940120897, 0.6502324818547323
1.1085896400378799, 0.6447212183644377
1.210053730887228, 0.6326359813834777
1.311538042434257, 0.623836607775757
1.4068752618861193, 0.6161347687033029
1.5144999252957558, 0.6051425727692359
1.6036765721873292, 0.5963476926498887
1.7020839675704944, 0.5875494424142615
1.7974279272549167, 0.5809428911328871
1.9019891549658103, 0.572142394153073
2.0034869469779606, 0.5655335961275118
2.1357067190008276, 0.5512465498434582
2.320234065809384, 0.5369404062338164
2.4278587292190204, 0.5259482102997495
2.5477775768189925, 0.5127609452951493
2.658499377090174, 0.505053489362228
2.7907056686479206, 0.4885758674960149
2.907561080549149, 0.47758030144566777
3.027466447684, 0.4622024608589079
3.1565890828453216, 0.44463067457370853
3.230360928219894, 0.43255554794158924
3.3010625976631616, 0.4215768324726431];

hold on

plot(regularDataAir(:,1),regularDataAir(:,2),'o','Color',string(lineCol(1)),'MarkerFaceColor',string(lineCol(1)),'MarkerEdgeColor','k','MarkerSize',6)
plot(minichannelData_Air(:,1),minichannelData_Air(:,2),'o','Color',string(lineCol(2)),'MarkerFaceColor',string(lineCol(2)),'MarkerEdgeColor','k','MarkerSize',6)
plot(minichannelData_o2(:,1),minichannelData_o2(:,2),'o','Color',string(lineCol(3)),'MarkerFaceColor',string(lineCol(3)),'MarkerEdgeColor','k','MarkerSize',6)

ylim([0 1])

% plot power density
yyaxis right
plot((i_flood./10000),(i_flood./10000).*V_flood,'--k','lineWidth',2);
hold on
plot((i_dry./10000),(i_dry./10000).*V_dry,'--r','lineWidth',2);
ylabel('Power Density (W cm^{-2})')

legend('Model Air - \delta_w = 299 \mum','Model Air - \delta_w = 0 \mum','Model O_2 - \delta_w = 0 \mum','Exp Air - Serpentine','Exp Air - Serpentine Micro-channel','Exp O_2 - Serpentine Micro-channel')

legend('Model Air - \delta_w = 209 \mum','Model Air - \delta_w = 0 \mum','Exp Air - Serpentine','Exp Air - Serpentine Micro-channel')
legend('\delta_w = 209 \mum','\delta_w = 0 \mum','Power Density \delta_w = 209 \mum','Power Density \delta_w = 0 \mum')

data_dry = [i_dry./10000 V_dry C_dry];
data_flood = [i_flood./10000 V_flood C_flood];

writematrix(data_dry,'dry209.csv')
writematrix(data_flood,'flood209.csv')

set(gca,'YColor','k')
set(gca,'color','none')

end
