%% Compare asymptotic estimates to numerics

% Here we estimate effective pressure and effective porosity, and then meltwater flux 
% from analytics presented in Ranganathan and others (submitted)
% and compare to results from numerics presented in Meyer et al. 2018

clear all;
close all;

% top-level parameters
theta = 1;
H = 200;
%H = 1000;
num_pts = 256;
%num_pts = 100;
z = linspace(0,H,num_pts);
z = z./H; % non dimensionalize

% thermomechanical parameters
Tm = 273.15;
Ts = 273.15-1;
dT = Tm-Ts;
K = 2.1; % W m^-1 K^-1
cp = 2009; % J kg^-1 K^-1

% rheological parameters
n = 3;
A = 2.4e-24; % Pa^-3 s^-1

% temperate-related parameters
kappa = 0.4416;
alpha = 2;
N0 = 1;
delta = 0.0023;
rho_water = 1000; % kg m^-3
rho_ice = 916; % kg m^-3;

% input parameters
Br = 22.4919;
%Br = 6;
Pe = 1.1115;

%(K.*dT./(rho_w.*L.*mean(Hvec)))

load('numerics_results.mat')

[T,zct,sr,crit_sr] = findIceTemperature(Br,Pe,theta,H,z,Tm,Ts,K,A,n,dT);
zct = zct./H; % non dimensionalize

[porosity_outer] = findOuterPorosityProfile(-Pe,Br,zct,z,kappa,alpha);
[effpressure_outer] = findOuterPressureProfile(Br,-Pe,porosity_outer,kappa,alpha);

[porosity_inner_0] = findInnerPorosity0Profile(porosity_outer);
[effpressure_inner] = findInnerPressureProfile(effpressure_outer,N0,alpha,kappa,-Pe,porosity_inner_0,z,zct,delta);
[porosity_inner_1] = findInnerPorosity1Profile(Br,-Pe,porosity_inner_0,z,zct,delta,effpressure_outer,alpha,kappa,N0);
[porosity_inner] = findInnerPorosityProfile(porosity_inner_0,porosity_inner_1,delta,z,zct);

[effpressure_composite] = findCompositePressureProfile(effpressure_outer,effpressure_inner);
[porosity_composite] = findCompositePorosityProfile(porosity_outer,porosity_inner,Br,kappa,alpha,-Pe,zct,z);

[flux] = findMeltwaterFlux(effpressure_composite,porosity_composite,alpha,kappa,z,delta);

% plotting
cmap = colorcet('l9','N',10,'reverse',0);

figure;
subplot(2,2,1)
%rectangle('Position',[-1 0 2 0.5],'FaceColor',[0.9 0.9 0.9])
hold on
scatter(Temperature(1:10:end),z(1:10:end),30,cmap(1,:),'d','filled')
plot((T-Tm)./dT,z,'LineWidth',2,'Color',cmap(5,:))
plot([-1 1],zct*ones(size([0 18])),'LineWidth',2,'Color',[0.7 0.7 0.7])
grid on
colormap(cmap)
xlabel('Temperature')
ylabel('Height')
xlim([-1 0])
ylim([0 1])
yticks([0 0.2 0.4 0.6 0.68 0.8 1])
set(gca,'FontSize',12,'FontWeight','b','GridColor','r');
legend('Numerics','Model (Meyer+Minchew 2018)')

EffP(z>=zct) = NaN;
effpressure_inner(z>=zct) = NaN;
effpressure_outer(z>=zct) = NaN;
effpressure_composite(z>=zct) = NaN;

subplot(2,2,2)
%rectangle('Position',[0 0 2 0.5],'FaceColor',[0.9 0.9 0.9])
scatter(EffP(1:10:end),z(1:10:end),30,cmap(1,:),'d','filled')
hold on
plot(effpressure_inner,z,'LineWidth',2,'Color',cmap(5,:))
plot(effpressure_outer,z,'--','LineWidth',2,'Color',cmap(5,:))
plot(effpressure_composite,z,'LineWidth',2,'Color',cmap(8,:))
plot([0 18],delta.^0.5*ones(size([0 18])),'LineWidth',2,'Color',[0.4 0.4 0.4])
grid on
colormap(cmap)
xlabel('Effective Pressure')
xlim([0 18])
xticks([0 2 4 6 8 10 12 14 16 18])
ylim([0 zct])
yticks([0 0.2 0.4 0.6 0.68])
set(gca,'FontSize',12,'FontWeight','b','GridColor','r');
legend('Numerics','Inner Solution','Outer Solution','Composite Solution')

Porosity(z>=zct) = NaN;
porosity_inner(z>=zct) = NaN;
porosity_outer(z>=zct) = NaN;
porosity_composite(z>=zct) = NaN;

subplot(2,2,3)
%rectangle('Position',[0 0 1 0.5],'FaceColor',[0.9 0.9 0.9])
scatter(Porosity(1:10:end),z(1:10:end),30,cmap(1,:),'d','filled')
hold on
plot(porosity_inner,z,'LineWidth',2,'Color',cmap(5,:))
plot(porosity_outer,z,'--','LineWidth',2,'Color',cmap(5,:))
plot(porosity_composite,z,'LineWidth',2,'Color',cmap(8,:))
plot([0 6],delta.^0.5*ones(size([0 6])),'LineWidth',2,'Color',[0.4 0.4 0.4])
grid on
colormap(cmap)
xlabel('Porosity (%)')
ylabel('Height')
xlim([0 6])
xticks([0 3 6])
ylim([0 zct])
yticks([0 0.2 0.4 0.6 0.68])
set(gca,'FontSize',12,'FontWeight','b','GridColor','r');
legend('Numerics','Inner Solution','Outer Solution','Composite Solution')

flux(z>=zct) = NaN;

subplot(2,2,4)
scatter(qdrain,0,100,cmap(1,:),'d','filled')
hold on
plot(flux,z,'LineWidth',2,'Color',cmap(5,:))
plot([-12 0],delta.^0.5*ones(size([-12 0])),'LineWidth',2,'Color',[0.4 0.4 0.4])
grid on
colormap(cmap)
xlabel('Flux')
xlim([-12 0])
xticks([-12 -10 -8 -6 -4 -2 0])
ylim([0 zct])
yticks([0 0.2 0.4 0.6 0.68])
set(gca,'FontSize',12,'FontWeight','b','GridColor','r');
legend('Numerics','Asymptotic Estimate')
