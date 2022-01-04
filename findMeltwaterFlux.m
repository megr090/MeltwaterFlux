function [flux] = findMeltwaterFlux(effpressure_composite,porosity_composite,alpha,kappa,z,delta)

dNdz = gradient(effpressure_composite)./gradient(z);
flux = (kappa.*porosity_composite.^alpha).*(-1+delta.*dNdz);


end