function [effpressure_outer] = findOuterPressureProfile(Br,Pe,porosity_outer,kappa,alpha)

effpressure_outer = (-Br.*kappa.*alpha.*porosity_outer.^(alpha-2))./(Pe-kappa.*alpha.*porosity_outer.^(alpha-1));

end