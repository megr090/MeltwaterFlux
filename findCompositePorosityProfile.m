function [porosity_composite] = findCompositePorosityProfile(porosity_outer,porosity_inner,Br,kappa,alpha,Pe,zct,z)

y1 = -Br.*zct./(kappa.*alpha.*porosity_outer(1).^(alpha-1)-Pe);

porosity_composite = porosity_outer + porosity_inner - (porosity_outer(1)+(z./zct).*y1);
porosity_composite(z>zct) = 0;


end

