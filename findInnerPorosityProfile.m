function [porosity_inner] = findInnerPorosityProfile(porosity_inner_0,porosity_inner_1,delta,z,zct)

porosity_inner = porosity_inner_0+porosity_inner_1.*delta.^(0.5);
porosity_inner(z>zct) = 0;

end

