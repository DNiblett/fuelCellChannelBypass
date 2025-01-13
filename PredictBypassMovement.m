%% Critical Air Velocity Required To Remove Liquid Water Plug in Serpentine By-pass Flow Field %%
% Collaboration between Newcastle University (Dr Daniel Niblett)
% and University of New South Wales (Dr Quentin Meyer)
% Code developed 24/09/2024

clc
clear
plotvariables = 1; % 1 for plotting figures

%% Predicting Critical Air Velocity for Bypass Size vs Distance from Serpentine Bend

for makingDesignSpace = 1
% Serpentine channel conditions
k_1 = 1e-8; % permeability of square channel (non-developed)

Hc = 0.001; % height of main channel
Wc = 0.001; % width of main channel
L = 0.01;  % length from the start to the end of the minichannel through the serpentine bend
cornerAngle = pi()/4; % half corner angle for serpentine bend
Lt = [1e-6:100e-6:0.02]; % Length from start to end of minichannel (i.e. Y-axis in figure 1) (m)
sigma = 0.072; % surface tension
theta_GDL = 120.*pi()./180; % contact angle gas diffusion layer
theta_w = 100.*pi()./180; % contact angle channel wall

% preallocate arrays
Hsave=[]; 
Usave =[];
Lsave = [];
Psave = [];

for t = 1:size(Lt,2)

% Minichannel conditions
H = 100e-6;
W = 100e-6;
L = Lt(t);
H = [50e-6:10e-6:300e-6]; % list of minichannel heights (I.e. X-axis in figure 1) (m)
W = H; % isotropic width to height

% Contact angle conditions
theta_a = 150*pi()/180;
theta_r = 85*pi()/180;

% Fluid properties air (to match that used in OpenFOAM simulations)
viscosity = 1e-5;
density = 1;
surfaceTension = 0.072;

% constants of critical velocity equation
a = 1.*density.*(1-cos((3/4).*(pi()-cornerAngle)));
bp=(Wc)/2;
cp=(Hc)/2;
betap=((1/3) - (64/(pi()^5)).*((cp./bp).*tanh((pi().*bp)./ (2.*cp))));
b = (viscosity.*L)./(k_1);
c = (-2.*sigma.*(sin(theta_w- pi()/2)+sin(theta_GDL - pi()./2))./(H));

% critical velocity using quadratic formula
U = (-b+sqrt(b.^2 - 4.*a.*c)) ./ (2.*a);

% estimate pressure drop
pressureDropCorner = 1.*U.^2 .*density.*(1-cos((3/4).*(pi()-cornerAngle)));
DPsinglephase = (U.*viscosity.*L)./(k_1);
DP = DPsinglephase + pressureDropCorner;

i
Hsave = [Hsave;H'];
Usave = [Usave;real(U')];
Lsave = [Lsave;[ones(size(U)).*L]'];
Psave = [Psave;DP'];

end

if plotvariables == 1
    
figure(1)


% plot pressure on the second y axis as line
yyaxis right
R = H./((sin(theta_w- (pi()/2))+sin(theta_GDL - pi()/2)));
Pc = 2*(sigma)./R;
plot(H.*1e+6,Pc,'-k','linewidth',2)
ylabel('Exit Capillary Pressure (Pa)')
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
hold on
line([100 100],[0 2000],'linestyle','--','color','k','linewidth',2)

yyaxis left

% convert U into L/min
Usave = Usave.*(0.001.^2).*60000;
DT = delaunay([Hsave,Lsave]);
trisurf(DT,Hsave.*1e+6,Lsave.*1e+3,Usave,'linestyle','none');
view(2)
shading interp
colormap('jet')
%cmap = cmocean('balance')
%cmap = cmocean('balance')
%colormap(cmap)
colorbar
set(gca,'fontsize',16)
xlabel('Lateral bypass width Size (\mum)')
ylabel('Distance From start to end of serpentine (mm)')
ylabel('Serpentine Detour Length (mm)')
%xlim([0.1 max(Hsave)*1e+3])
hcb=colorbar;
hcb.Title.String = "Critical Velocity (m/s)";
hcb.Title.String = "Critical Flow Rate (L min^{-1})";
xlim([min(Hsave).*1e+6 max(Hsave).*1e+6])
box off
grid off

% plot pressure on the second y axis as line
yyaxis right
R = H./((sin(theta_w- (pi()/2))+sin(theta_GDL - pi()/2)));
Pc = 2*(sigma)./R
plot(H.*1e+3,Pc,'-k','linewidth',2)
ylabel('Exit Capillary Pressure (Pa)')
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
hold on
%line([100 100],[0 200],'linestyle','--','color','k','linewidth',2)


%Also plot new dataset ontop with markers
data = [16.9 20 1 1 120 70;
        12.9 20 1 1 120 70;
        8.9 20 1 0 120 70;
        5.0 20 1 0 120 70;
        16.9 17 1 0 120 70;
        12.9 17 1 0 120 70;
        8.9 17 0 0 120 70;
        5.0 17 0 0 120 70;
        16.9 15 0 0 120 70;
        12.9 15 0 0 120 70;
        8.9 15 0 0 120 70;
        5.0 15 0 0 120 70;
        16.9 10 0 0 120 70;
        12.9 10 0 0 120 70;
        8.9 10 0 0 120 70;
        5.0 10 0 0 120 70;
        16.9 15 1 0 150 85;
        12.9 15 1 0 150 85;
        8.9 15 0 0 150 85;
        5.0 15 0 0 150 85;];

% data from the CFD for distance vs detach or not.
data = [16.9 20 0 0 120 100;
        12.9 20 0 0 120 100;
        8.9 20 0 0 120 100;
        5.0 20 0 0 120 100];

set(gca,'color','none')
end


%% Or instead plot the channel distance vs velocity and then make line for critical velocity

% Main Channel Height and Width
Hc = 0.001;
Wc = 0.001;

% Minichannel Height & Width
H = 100e-6;
W = 100e-6;

L = [0.001:0.001:0.02];

% Contact angle conditions
theta_a = 120*pi()/180;
theta_r = 70*pi()/180;

%theta_a = 150*pi()/180;
%theta_r = 85*pi()/180;

% Fluid properties air
viscosity = 1e-5;
density = 1;
surfaceTension = 0.072;

% constants of critical velocity equation
a = 1.*density.*(1-cos((3/4).*(pi()-cornerAngle)));
bp=(Wc)/2;
cp=(Hc)/2;
betap=((1/3) - (64/(pi()^5)).*((cp./bp).*tanh((pi().*bp)./ (2.*cp))));
b = (viscosity.*L)./(k_1);
c = (-2.*sigma.*(sin(theta_w- pi()/2)+sin(theta_GDL - pi()./2))./(H));
U = (-b+sqrt(b.^2 - 4.*a.*c)) ./ (2.*a);
pressureDropViscous = U.*(viscosity.*L)./((Hc./2).*( (1/3) - ( 64.*Hc./(Wc.*pi().^5) ).*tanh(pi().*Wc./(2.*Hc))));
pressureDropCorner = 1.*U.^2 .*density.*(1-cos((3/4).*(pi()-cornerAngle)));
DPsinglephase = (U.*viscosity.*L)./(k_1);
DP = DPsinglephase + pressureDropCorner;


% 
% % plot L vs U
% figure(2)
% plot(L.*1e+3,U,'k','linewidth',2)
% xlabel('Distance from start to end of mini-channel (mm)')
% ylabel('Critical Velocity (m/s)')
% set(gca,'fontsize',16)
% hold on

% plot L vs U
figure(2)
plot(L.*1e+3,U.*(0.001^2).*60000,'k','linewidth',2)
xlabel('Distance from start to end of mini-channel (mm)')
ylabel('Critical Flow Rate (L min^{-1})')
set(gca,'fontsize',16)
hold on

data = [0.005 1 1.5;
0.009 2 1.5;
0.013 2 1.5;
0.017 2 1.5;
0.005 1 1.2;
0.009 1 1.2;
0.013 1 1.2;
0.017 1 1.2;
0.005 1 0.9;
0.009 1 0.9;
0.013 1 0.9;
0.017 1 0.9;
0.005 2 1.8;
0.009 2 1.8;
0.013 2 1.8;
0.017 2 1.8;];

xl = data(data(:,2)==2);
yl = (xl./xl).*data(data(:,2)==2,3);
hold on
plot(xl.*1e+3,yl,'ko','MarkerFaceColor','g','MarkerSize',10)

xl = data(data(:,2)==1);
yl = (xl./xl).*data(data(:,2)==1,3);
hold on
plot(xl.*1e+3,yl,'kv','MarkerFaceColor','r','MarkerSize',10)


%idv = (data(:,4) == 0 & data(:,5) == 120)
%idv = (data(:,3) == 0 & data(:,5) == 150)

%xl = data(idv,1);
%yl = data(idv,2);
%hold on
%plot(xl,yl,'kv','MarkerFaceColor','r','MarkerSize',10)


%legend('Critical Velocity','CFD detach','CFD no-detach')

end

%% Test prediction of pressure drop along channel vs OpenFOAM CFD %%

for comparePressureDropPredictions = 1

% Main Channel Height and Width
Hc = 0.001;
Wc = 0.001;

% Minichannel Height & Width
H = 100e-6;
W = 100e-6;

L = [0.001:0.001:0.02];

U = 25;

% constants of critical velocity equation
a = 1.5.*density.*(1-cos((3/4).*(pi()-cornerAngle)));
bp=(Wc)/2;
cp=(Hc)/2;
betap=((1/3) - (64/(pi()^5)).*((cp./bp).*tanh((pi().*bp)./ (2.*cp))));
%b = (viscosity.*L)./(cp.^2 .* betap);
%c = (-2.*sigma.*(sin(theta_w- pi()/2)+sin(theta_GDL - pi()./2))./(H));
%U = (-b+sqrt(b.^2 - 4.*a.*c)) ./ (2.*a);
pressureDropViscous = U.*(viscosity.*L)./((Hc./2).*( (1/3) - ( 64.*Hc./(Wc.*pi().^5) ).*tanh(pi().*Wc./(2.*Hc))));
pressureDropCorner = 1.5.*U.^2 .*density.*(1-cos((3/4).*(pi()-cornerAngle)));
DPsinglephase = (U.*viscosity.*L)./(cp.^2 .* betap);
DP = DPsinglephase + pressureDropCorner;

first_channel_L = [0:0.001:0.008];
L = first_channel_L;
pressureDropCorner = 1.5.*U.^2 .*density.*(1-cos((3/4).*(pi()-cornerAngle)));
DP1 = (U.*viscosity.*L)./(cp.^2 .* betap);

second_channel_L = [0:0.001:0.002];
L = second_channel_L;
DP2 = (U.*viscosity.*L)./(cp.^2 .* betap);
pressureDropCorner = 1.5.*U.^2 .*density.*(1-cos((3/4).*(pi()-cornerAngle)));

third_channel_L = [0:0.001:0.008];
L = third_channel_L;
pressureDropCorner = 1.5.*U.^2 .*density.*(1-cos((3/4).*(pi()-cornerAngle)));
DP3 = (U.*viscosity.*L)./(cp.^2 .* betap);

totalL = [first_channel_L';[second_channel_L+first_channel_L(end)]';[third_channel_L+second_channel_L(end)+first_channel_L(end)]'];
P1 = 1250-DP1;
P2 = P1(end)-DP2-pressureDropCorner./2;
P3 = P2(end)-(pressureDropCorner./2)-DP3;

totalP = [P1';P2';P3'];

% load data
data = csvread('./pressure_channel_data/pressure_channel.csv',1);

% plot L vs U
figure(2)
plot(totalL.*1e+3,totalP,'k','linewidth',2)
xlabel('Distance from start of main channel to end (mm)')
ylabel('Pressure (Pa)')
set(gca,'fontsize',16)
hold on
plot(data(:,9).*1e+3,data(:,10),'b','linewidth',2)
legend('Analytical','StreamLine CFD Extraction')

end


%% Test prediction of pressure drop along channel vs OpenFOAM CFD %%
% This time with corrected square channel permeability

for comparePressureDropPredictionCorrected = 1
k_2 = k_1;
k_1 = 1e-8;
% Main Channel Height and Width
Hc = 0.001;
Wc = 0.001;

% Minichannel Height & Width
H = 100e-6;
W = 100e-6;
C_corner = 1;
cornerAngle = pi()/4; % half corner angle for serpentine bend

L = [0.001:0.001:0.02];
U = 25;

% constants of critical velocity equation
a = C_corner.*density.*(1-cos((3/4).*(pi()-cornerAngle)));
bp=(Wc)/2;
cp=(Hc)/2;
betap=((1/3) - (64/(pi()^5)).*((cp./bp).*tanh((pi().*bp)./ (2.*cp))));
%b = (viscosity.*L)./(cp.^2 .* betap);
%c = (-2.*sigma.*(sin(theta_w- pi()/2)+sin(theta_GDL - pi()./2))./(H));
%U = (-b+sqrt(b.^2 - 4.*a.*c)) ./ (2.*a);
pressureDropViscous = U.*(viscosity).*L;
pressureDropCorner = C_corner.*U.^2 .*density.*(1-cos((3/4).*(pi()-cornerAngle)));
DPsinglephase = (U.*viscosity.*L)./(cp.^2 .* betap);
DP = DPsinglephase + pressureDropCorner;

first_channel_L = [0:0.001:0.008];
L = first_channel_L;
pressureDropCorner = C_corner.*U.^2 .*density.*(1-cos((3/4).*(pi()-cornerAngle)));
DP1 = (U.*viscosity.*L)./k_1;

second_channel_L = [0.001:0.001:0.002];
L = second_channel_L;
pressureDropCorner = C_corner.*U.^2 .*density.*(1-cos((3/4).*(pi()-cornerAngle)));
DP2 = (U.*viscosity.*L)./k_1;
DP2 = DP2.*0;
DP2 = pressureDropCorner./size(L,2);
DP2 = [ones(size(L,2),1).*DP2]';
DP2 = cumsum(DP2);
pressureDropCorner = C_corner.*U.^2 .*density.*(1-cos((3/4).*(pi()-cornerAngle)));

third_channel_L = [0.001:0.001:0.008];
L = third_channel_L;
pressureDropCorner = C_corner.*U.^2 .*density.*(1-cos((3/4).*(pi()-cornerAngle)));
DP3 = (U.*viscosity.*L)./((k_1+k_2)./2);

totalL = [first_channel_L';[second_channel_L+first_channel_L(end)]';[third_channel_L+second_channel_L(end)+first_channel_L(end)]'];
P1 = 1235-DP1;
P2 = P1(end)-DP2;
P3 = P2(end)-DP3;

totalP = [P1';P2';P3'];

% load data
data = csvread('./pressure_channel_data/pressure_channel.csv',1);

% plot L vs U
figure(2)
plot(totalL.*1e+3,totalP,'k','linewidth',2)
xlabel('Distance from start of main channel to end (mm)')
ylabel('Pressure (Pa)')
set(gca,'fontsize',16)
hold on
plot(data(:,9).*1e+3,data(:,10),'b','linewidth',2)
%plot(data2(:,9).*1e+3,data2(:,10),'b','linewidth',2)
legend('Analytical','StreamLine CFD Extraction')


end


%% Analytical Prediction of Single Serpentine Bend Pressure Drop vs air flow rate along with critical velocity%

for predictPressureVsFlowRate = 1
u = [0.1:0.1:30];
L = 0.016;
cornerAngle = pi()/(4);
density = 1.0
N = 1.5;
for t = 1:size(u,2);
U = u(t);
pressureDropViscous = U.*(viscosity.*L)./((Hc./2).*( (1/3) - ( 64.*Hc./(Wc.*pi().^5) ).*tanh(pi().*Wc./(2.*Hc))));
%pressureDropCorner = 2.*U.^2 .*density.*(1-cos((3/4).*(pi()-cornerAngle)));
pressureDropCorner = N.*U.^2 .*density.*(1-cos((3/4).*(pi()-cornerAngle)));
DPsinglephase = (U.*viscosity.*L)./(cp.^2 .* betap);
%DPsinglephase2 = (2.*U.*viscosity.*0.002)./(cp.^2 .* betap);
DP = DPsinglephase + pressureDropCorner;
DPsave(t,1) = DP;
end

u2flow = @(u) u.*(0.001.^2).*6e+7./1000;
Q_a = u2flow(u);
figure(3)
plot(Q_a,DPsave,'-k','linewidth',2)
xlabel('Air Flow Rate (L min^{-1})')
ylabel('Pressure Drop (Pa)')

data = [1 9.622;5 75;10 240.213;15 458.18;20 748.56] % from CFD data
data2 = [25 1173;30 1679] % from CFD data
hold on
plot(u2flow(data(:,1)),data(:,2),'vk','MarkerSize',10,'MarkerFaceColor','r')
plot(u2flow(data2(:,1)),data2(:,2),'ok','MarkerSize',10,'MarkerFaceColor','g')
set(gca,'fontsize',16)
grid off
line([u2flow(23.6) u2flow(23.6)],[0 1800],'color','g','linestyle','--','linewidth',2)
%line([22 22],[0 1800],'color','k','linestyle','--','linewidth',2)

%line([25 25],[0 1800],'color','g','linestyle','--','linewidth',2)
%line([30 30],[0 1800],'color','g','linestyle','--','linewidth',2)
legend('\DeltaP Analytical Prediction','\DeltaP CFD Extraction - No Plug Removal','\DeltaP CFD Extraction - Plug Removal','Critical Q_a Prediction')

ylim([0 1800])

end


%% Breakthrough Capillary Pressure Required to exit the bypass channel %%

for predictBreakThroughPressure = 1
H = 100e-6;
theta_w = [90:1:180].*pi()/180;
theta_s = [90:1:180].*pi()/180;
[X,Y] = meshgrid(theta_w,theta_s);
sigma = 0.072;
R = H./(sin(X - pi()/2)+sin(Y - pi()/2));
Pc = 2.*sigma./R;

figure(4)
%plot(Pc)
surf(X.*180./pi(),Y.*180./pi(),Pc,'linestyle','none')
view(2)
h = colorbar;
%colormap('jet')
cmap = cmocean('matter');
colormap(cmap)

h.Label.String = 'Breakthrough Pressure (Pa)';
hold on
view(2)
ylabel('GDL Contact Angle ({\circ})')
xlabel('Channel Contact Angle ({\circ})')
set(gca,'fontsize',16)
xlim([min(min(X)).*180./pi() 180])
ylim([min(min(X)).*180./pi() 180])

% also plot the entry breakthrough pressure for our system
H = [50e-6:1e-6:1000e-6];
theta_w = 100.*pi()/180;
theta_s = 120.*pi()/180;

sigma = 0.072;
R = H./(sin(theta_w - pi()/2)+sin(theta_s - pi()/2));
Pc = 2.*sigma./R;

figure(5)
plot(H.*1e+6,Pc,'k')
set(gca,'fontsize',16)
ylabel('Breakthrough Pressure (Pa)')
xlabel('Size of mini-channel (\mum)')
hold on

H = 100e-6;
theta_w = [0:1:180].*pi()/180;
theta_s = [0:1:180].*pi()/180;
[X,Y] = meshgrid(theta_w,theta_s);
sigma = 0.072;
R = H./(sin(X)+sin(Y - pi()/2));
Pc = 2.*sigma./R;

figure(6)

surf(X.*180./pi(),Y.*180./pi(),Pc,'linestyle','none')
view(2)
h = colorbar;
colormap('jet')
h.Label.String = 'Breakthrough Pressure (Pa)'
hold on
view(2)
ylabel('GDL Contact Angle')
xlabel('Channel Contact Angle')
set(gca,'fontsize',16)
xlim([0 180])
ylim([0 180])

H = [50e-6:1e-6:1000e-6];
theta_w = 100.*pi()/180;
theta_s = 120.*pi()/180;

sigma = 0.072;
R = H./((sin(theta_w- (pi()/2))+sin(theta_s - pi()/2)));
Pc = 2*(sigma)./R


figure(7)
plot(H.*1e+6,Pc,'linewidth',2)
set(gca,'fontsize',16)
ylabel('Breakthrough Pressure (Pa)')
xlabel('Size of mini-channel (\mum)')
set(gca,'fontsize',16) 
hold on

end


