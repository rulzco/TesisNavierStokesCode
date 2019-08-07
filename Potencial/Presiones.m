function [Cp,p]=Presiones(u,v,Vinf,dinf,gamma,pinf,p0,d0,H0)
%----------------------------VARIABLES--------------------------------%
%p= Presión
%Cp= Coeficiente de resistencia
%---------------------------------------------------------------------%
%Aplicamos la fórmula (3.5).
    d=d0.*(1-(u.^2+v.^2)/(2*H0)).^(1/(gamma-1));
    p=p0*(d/d0).^gamma;
    Cp=real(2*(p-pinf)./(dinf*Vinf^2));
end