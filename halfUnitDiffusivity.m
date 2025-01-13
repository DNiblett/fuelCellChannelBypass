%% Fuel Cell Cathode - Half Channel, Half Rib Effective Diffusivity Model %%
% Equation developed analytically by Dr Daniel Niblett 20/11/2023

%.............
%      :     |
% CHNL :     |
%______:  G  | C
%      |  D  | L
%  RIB |  L  |
% .....|.....|

Wc = 500e-6; % Half Channel Width
Wr = 500e-6; % Half Rib Width
D = 2e-5;
porosity = 0.78;
tau = 1.5;
delta_g = 209e-6;
delta_w = [0e-6:1e-6:209e-6];

Dratio = delta_g./((Wc+Wr).* ((delta_w./Wc) + (sqrt((Wc/2).^2 + (delta_g -delta_w).^2)./sqrt((Wc+Wr).^2 + (delta_g -delta_w).^2)    )   ) );
Dporous = D.*porosity.^tau
Dactual = Dporous.*Dratio;
% 
figure(1)
plot((delta_w.*1e+6),Dactual./Dporous,'k','LineWidth',2)
xlabel('Water film thickness (\mum)')
ylabel('Effective Diffusivity')
set(gca,'fontsize',16)
xlim([0 max(delta_w).*1e+6])

% figure(2)
% plot((delta_w)./(delta_g),Dactual./Dporous,'k','LineWidth',2)
% xlabel('Water film thickness (\mum)')
% ylabel('Effective Diffusivity')
% set(gca,'fontsize',16)

data = [[delta_w]'.*1e+6 [Dactual./Dporous]'];

writematrix(data,'waterThickness_Diffusivity.csv')



