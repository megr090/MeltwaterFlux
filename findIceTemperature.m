function [T,zct,sr,crit_sr] = findIceTemperature(Br,Pe,theta,H,z,Tm,Ts,K,A,n,dT)

z = z.*H;

% solve for strain rate
sr = ((K.*dT.*Br)./(2.*H.^2.*theta)).^(n/(n+1)).*A.^(1./(n+1));

% find critical strain rate
crit_sr = ((0.5.*Pe.^2)./(Pe-1+exp(-Pe))).^(n./(n+1)).*((K.*dT)./(A^(-1/n).*H.^2.*theta)).^(n/(n+1));

% find the thickness of the temperate zone
if sr > crit_sr
    zetaH = 1-(Pe/(Br))-(1/Pe).*(1+real(lambertw(-exp((-Pe.^2)./(Br)-1))));
else
    zetaH = 0;
end

zct = zetaH.*H;

% find the temperature profile
for i=1:length(z)
    if and(z(i)>=zetaH.*H,z(i)<=H)
        T(i) = Ts + dT*(Br./Pe).*(1-(z(i)./H)+(1/Pe).*exp(Pe.*(zetaH-1))-(1/Pe).*exp(Pe.*(zetaH.*H-z(i))./H));
    else
        T(i) = Tm;
    end
end

end

