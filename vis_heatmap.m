close all
clear all
mu_f = 5e-5;
delta_rhog = 1000*9.81;
eta_phi = 5e15;
phi_0 = 0.02;
ly = 70;
lx = ly/2;
ny = 256*2-1;
nx = 256-1;
dy = ly/(ny) ;
k_s = 5e-17;
delta_c = sqrt(k_s*eta_phi/mu_f);
perm = 5e-15; % m^2
mu_f = 5e-5; % viscosity of fluid, Pas
eta_phi = 5e15; % Bulk viscosity, Pas
delta_rhog = 1000*9.81; % Buoyancy phase differential
tau_c = ((eta_phi*mu_f/perm)^0.5)/delta_rhog; % Characteristic time
FS = 24;
path={'/Users/kmagno/Documents/porosity_update/np30/'};
import_path = string(path{1});
A = dir(strcat(import_path,'qDys*.csv'));
names = {A(:).name};
B=natsort(names);
Z = zeros([512 255]);
for jj = 1:length(B)
    T = readmatrix(strcat(import_path,string(B{jj})));
    T = rot90(T);
    figure(3),set(gcf,'Position',[45 99 1066 663]);
    imagesc(T);colorbar; crameri('batlowW'); colorbar;crameri('batlowW'); hold on;
    set(gca, 'clim', [0 200]);
    ax2 = gca;
    ax2.FontSize = FS;
    ax2.TickLabelInterpreter = 'latex';
end
