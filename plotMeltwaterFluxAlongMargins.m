function [done] = plotMeltwaterFluxAlongMargins(x,flux,pressure_composite,porosity_composite,N_dim,flux_dim,N_dimcomp,SRvec,Hvec,SMBvec,Tsvec,epsilon,N0,zct,flux_dimcomp)

flux_dim = flux_dim.*(240^2);

% compute rate of shear heating
n = 3;
A = 2.4e-24; % Pa s^-1
W = 2.*A.^(-1/n).*SRvec.^((n+1)./n);
z = linspace(0,1,100);

figure;
subplot(4,1,1)
semilogy(x,SRvec,'LineWidth',2,'Color','k')
set(gca,'FontSize',12,'FontWeight','b','GridColor','r');
ylabel('Strain Rate')

subplot(4,1,2)
plot(x,Hvec,'LineWidth',2,'Color','k')
set(gca,'FontSize',12,'FontWeight','b','GridColor','r');
ylabel('H')

subplot(4,1,3)
plot(x,SMBvec,'LineWidth',2,'Color','k')
set(gca,'FontSize',12,'FontWeight','b','GridColor','r');
ylabel('SMB')

subplot(4,1,4)
plot(x,Tsvec,'LineWidth',2,'Color','k')
set(gca,'FontSize',12,'FontWeight','b','GridColor','r');
xlabel('Distance Upstream (km)')
ylabel('Ts')

figure;
subplot(3,1,1)
imagesc(x./1e3,z(end:-1:1),real(pressure_composite))
hold on;
plot(x./1e3,1-zct,'LineWidth',2,'Color','k')
c = colorbar;
colormap(colorcet('l17'));
set(get(c,'label'),'string','N (-)')
set(gca,'FontSize',12,'FontWeight','b','GridColor','r');
ylabel('z/H')

subplot(3,1,2)
imagesc(x./1e3,z(end:-1:1),real(porosity_composite))
hold on;
plot(x./1e3,1-zct,'LineWidth',2,'Color','k')
c = colorbar;
colormap(colorcet('l17'));
set(get(c,'label'),'string','Porosity')
set(gca,'FontSize',12,'FontWeight','b','GridColor','r');
ylabel('z/H')

subplot(3,1,3)
imagesc(x./1e3,z(end:-1:1),real(flux))
hold on;
plot(x./1e3,1-zct,'LineWidth',2,'Color','k')
c = colorbar;
colormap(colorcet('l17'));
set(get(c,'label'),'string','Flux (-)')
xlabel('Distance Upstream (km)')
ylabel('z/H')
set(gca,'FontSize',12,'FontWeight','b','GridColor','r');

figure;
subplot(3,1,1)
imagesc(x./1e3,z(end:-1:1),real(N_dim./1e3))
hold on;
plot(x./1e3,1-zct,'LineWidth',2,'Color','k')
c = colorbar;
colormap(colorcet('l17'));
set(get(c,'label'),'string','N (kPa)')
set(gca,'FontSize',12,'FontWeight','b','GridColor','r');
ylabel('z/H')

subplot(3,1,2)
imagesc(x./1e3,z(end:-1:1),real(porosity_composite))
hold on;
plot(x./1e3,1-zct,'LineWidth',2,'Color','k')
c = colorbar;
colormap(colorcet('l17'));
set(get(c,'label'),'string','Porosity (%)')
set(gca,'FontSize',12,'FontWeight','b','GridColor','r');
ylabel('z/H')

subplot(3,1,3)
imagesc(x./1e3,z(end:-1:1),real(flux_dim.*3.15e7))
hold on;
plot(x./1e3,1-zct,'LineWidth',2,'Color','k')
c = colorbar;
colormap(colorcet('l17'));
set(get(c,'label'),'string','Flux (m^3/yr)')
xlabel('Distance Upstream (km)')
ylabel('z/H')
set(gca,'FontSize',12,'FontWeight','b','GridColor','r');

done = 1;
end

