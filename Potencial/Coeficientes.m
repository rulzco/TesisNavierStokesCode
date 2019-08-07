function [L,D]=Coeficientes(N,M,Cp,X,Y,alf,cN)
% Inicializamos dos vectores en los cuales se pretende almacenar
% los coeficientes que se encuentran dentro de la integral de
% las ecuaciones de (3.15)
    ky=zeros(N,1);
    kx=zeros(N,1);
    xi=1:N;
    alfa=alf*pi/180;
    [g11,g22,g21,g12,J,dxdxi,dxdeta,dydxi,dydeta]=Tensor(X,Y,N,M);
    for i=1:N
        ky(i,1)=Cp(i,M)*dxdxi(i,M)/cN;
        kx(i,1)=Cp(i,M)*dydxi(i,M)/cN;
    end
% Integramos según la variable xi.
    CLy=trapz(xi,ky);
    CLx=-trapz(xi,kx);
% Descomponemos la fuerza en horizontal o vertical con el ángulo
% de ataque.
    L=-CLx*sin(alfa)+CLy*cos(alfa)
    D=CLx*cos(alfa)+CLy*sin(alfa)
end