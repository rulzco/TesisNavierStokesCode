function [g11,g22,g21,g12,J,dxdxi,dxdeta,dydxi,dydeta,A,B,C1]=Tensor(X,Y,N,M)
%----------------------------VARIABLES--------------------------------%
%J= Jacobiano
%g11,g12,g21,g22= Tensor métrico transformación directa
%g11I,g12I,g21I,g22I= Tensor métrico transformación indirecta
%dxdxi= Variación de la coordenada x con respecto a xi
%dxdeta= Variación de la coordenada x con respecto a eta
%dydxi= Variación de la coordenada y con respecto a xi
%dydeta= Variación de la coordenada y con respecto a eta
%---------------------------------------------------------------------%
% Para el cálculo de las variaciones de las coordenadas x e y
% se utilizan las discretizaciones de la tabla 4.3.1
    for i=1:N-1
        if i==1
            for j=2:M-1
                dxdxi(i,j)=(X(i+1,j)-X(N-1,j))/2;
                dxdeta(i,j)=(X(i,j+1)-X(i,j-1))/2;
                dydxi(i,j)=(Y(i+1,j)-Y(N-1,j))/2;
                dydeta(i,j)=(Y(i,j+1)-Y(i,j-1))/2;
            end
        else
            for j=2:M-1
                dxdxi(i,j)=(X(i+1,j)-X(i-1,j))/2;
                dxdeta(i,j)=(X(i,j+1)-X(i,j-1))/2;
                dydxi(i,j)=(Y(i+1,j)-Y(i-1,j))/2;
                dydeta(i,j)=(Y(i,j+1)-Y(i,j-1))/2;
             end
                dxdxi(i,M)=(X(i+1,M)-X(i-1,M))/2;
                dxdeta(i,M)=(X(i,M-2)-4*X(i,M-1)+3*X(i,M))/2;
                dydxi(i,M)=(Y(i+1,M)-Y(i-1,M))/2;
                dydeta(i,M)=(Y(i,M-2)-4*Y(i,M-1)+3*Y(i,M))/2;
                dxdxi(i,1)=(X(i+1,1)-X(i-1,1))/2;
                dxdeta(i,1)=(-3*X(i,1)+4*X(i,2)-X(i,3))/2;
                dydxi(i,1)=(Y(i+1,1)-Y(i-1,1))/2;
                dydeta(i,1)=(-3*Y(i,1)+4*Y(i,2)-Y(i,3))/2;
        end
    end
dxdxi(1,M)=(X(2,M)-X(N-1,M))/2;
dxdeta(1,M)=(X(1,M-2)-4*X(1,M-1)+3*X(1,M))/2;
dydxi(1,M)=(Y(2,M)-Y(N-1,M))/2;
dydeta(1,M)=(Y(1,M-2)-4*Y(1,M-1)+3*Y(1,M))/2;
dxdxi(1,1)=(X(2,1)-X(N-1,1))/2;
dxdeta(1,1)=(-X(1,3)+4*X(1,2)-3*X(1,1))/2;
dydxi(1,1)=(Y(2,1)-Y(N-1,1))/2;
dydeta(1,1)=(-Y(1,3)+4*Y(1,2)-3*Y(1,1))/2;
dydeta(N,:)=dydeta(1,:);
dydxi(N,:)=dydxi(1,:);
dxdeta(N,:)=dxdeta(1,:);
dxdxi(N,:)=dxdxi(1,:);
% Las fórmulas de los tensores métricos para las transformaciones
% directas o indirectas se encuentran en (2.21), (2.22) y (2.23).
J=dxdxi.*dydeta-dxdeta.*dydxi;
g11I=dxdxi.^2+dydxi.^2;
g12I=dxdxi.*dxdeta+dydxi.*dydeta;
g22I=dxdeta.^2+dydeta.^2;
g11=g22I./(J.^2);
g12=-g12I./(J.^2);
g22=g11I./(J.^2);
g21=g12;

C1=g11I;
A=g22I;
B=g12I;
end