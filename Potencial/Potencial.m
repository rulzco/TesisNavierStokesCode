function [phi,C,theta,IMA] = Potencial(d0,H0,gamma,Machinf,X,Y,Vinf,alf,N,M)
%----------------------------VARIABLES--------------------------------%
%d= Densidad
% Las variables ya definidas seguidas de una V o una H corresponden a:
% H: iteración en el nodo i '+' o '-' 1/2
% V: iteración en el nodo j '+' o '-' 1/2
%---------------------------------------------------------------------%
%----------------PARÁMETROS UTILIZADOS EN EL BUCLE--------------------%
[g11,g22,g21,g12,J,dxdxi,dxdeta,dydxi,dydeta]=Tensor(X,Y,N,M);
% Calculamos los coeficientes en los nodos intercalados.
% Fórmulas (4.14)
dxdxiV=zeros(N,M); dxdetaV=zeros(N,M);
dydxiV=zeros(N,M);dydetaV=zeros(N,M);
dxdxiH=zeros(N,M); dxdetaH=zeros(N,M);
dydxiH=zeros(N,M); dydetaH=zeros(N,M);
for i=1:N-1
for j=1:M-1
g11V(i,j)=0.5*(g11(i,j)+g11(i,j+1));
g12V(i,j)=0.5*(g12(i,j)+g12(i,j+1));
g22V(i,j)=0.5*(g22(i,j)+g22(i,j+1));
JV(i,j)=0.5*(J(i,j)+J(i,j+1));
dxdxiV(i,j)=0.5*(dxdxi(i,j)+dxdxi(i,j+1));
dxdetaV(i,j)=0.5*(dxdeta(i,j)+dxdeta(i,j+1));
dydxiV(i,j)=0.5*(dydxi(i,j)+dydxi(i,j+1));
dydetaV(i,j)=0.5*(dydeta(i,j)+dydeta(i,j+1));
end
end
for i=1:N-1
for j=1:M
g11H(i,j)=0.5*(g11(i,j)+g11(i+1,j));
g12H(i,j)=0.5*(g12(i,j)+g12(i+1,j));
g22H(i,j)=0.5*(g22(i,j)+g22(i+1,j));
JH(i,j)=0.5*(J(i,j)+J(i+1,j));
dxdxiH(i,j)=0.5*(dxdxi(i,j)+dxdxi(i+1,j));
dxdetaH(i,j)=0.5*(dxdeta(i,j)+dxdeta(i+1,j));
dydxiH(i,j)=0.5*(dydxi(i,j)+dydxi(i+1,j));
dydetaH(i,j)=0.5*(dydeta(i,j)+dydeta(i+1,j));
end
end
g21V=g12V;
g21H=g12H;
% Calculamos el ángulo theta de cada nodo. Como el arco tangente
% nos devuelve el valor del ángulo en los cuadrantes I y IV lo
% recalculamos según el siguiente código.
theta=atan(Y./X);
alfa=alf*pi/180;
for i=1:N
for j=1:M
if X(i,j)<=0 && Y(i,j)>=0
    theta(i,j)=pi-abs(theta(i,j));
elseif X(i,j)<=0 && Y(i,j)<0
theta(i,j)=abs(theta(i,j))+pi;
elseif X(i,j)>0 && Y(i,j)<0
theta(i,j)=2*pi-abs(theta(i,j));
end
end
end
theta(N,:)=2*pi;
theta(1,:)=0;
%----------------------------VALOR INICIAL----------------------------%
C=0.5;
phi=zeros(N,M);
UH=zeros(N,M);
VH=zeros(N,M);
UV=zeros(N,M);
VV=zeros(N,M);
ite=0; ddd=1; itemax=20000; tol=1.e-9;w=0.5;
while ddd>tol && ite<itemax
ite=ite+1;
phiold=phi;
%--------------------------FRONTERA EXTERIOR--------------------------%
% Para aplicar la fórmula (2.30) primero determinamos el arco tangente
% y lo distribuimos igual que en el caso de theta.
arcotan(:,1)=atan(sqrt(1-Machinf^2)*tan(theta(:,1)-alfa));
arcosen(:,1)=asin(sqrt(1-Machinf^2)*sin(theta(:,1)-alfa));
for i=1:N
if arcotan(i,1)>0&&arcosen(i,1)<0
arcotan(i,1)=arcotan(i,1)+pi;
elseif arcotan(i,1)<0&&arcosen(i,1)>0
arcotan(i,1)=pi- abs(arcotan(i,1));
elseif arcotan(i,1)<0&&arcosen(i,1)<0
if (theta(i,1)-alfa)>0
arcotan(i,1)=2*pi+arcotan(i,1);
end
end
end
phi(:,1)=Vinf*(X(:,1)*cos(alfa)+Y(:,1)*sin(alfa))+C*...
arcotan(:,1)/(2*pi);
%---------------------NODOS INTERNOS DE LA MALLA----------------------%
% Desarrollamos los parámetros en los nodos intercalados desarrollando
% las fórmulas (4.9) y (4.10).
for i=1:N-1
for j=1:M-1
if i==1&&j==M-1
PV(i,j)=0.25*(phi(i+1,j)-phi(N-1,j)+phi(i+1,j-1)...
-phi(N-1,j-1)+2*C);
elseif i==1&&j~=M-1
PV(i,j)=0.25*(phi(i+1,j+2)-phi(N-1,j+2)+phi(i+1,...
j+1)-phi(N-1,j+1)+2*C);
elseif i~=1&&j==M-1
PV(i,j)=0.25*(phi(i+1,j)-phi(i-1,j)+phi(i+1,j-1)...
-phi(i-1,j-1));
else
PV(i,j)=0.25*(phi(i+1,j+2)-phi(i-1,j+2)+phi(i+1,...
j+1)-phi(i-1,j+1));
end
if j==M-1
UV(i,j)=g11V(i,j)*PV(i,j)+g12V(i,j)*(phi(i,j)-...
phi(i,j-1));
VV(i,j)=g21V(i,j)*PV(i,j)+g22V(i,j)*(phi(i,j)-...
phi(i,j-1));
else
UV(i,j)=g11V(i,j)*PV(i,j)+g12V(i,j)*(phi(i,j+2)-...
phi(i,j+1));
VV(i,j)=g21V(i,j)*PV(i,j)+g22V(i,j)*(phi(i,j+2)-...
phi(i,j+1));
end
end
end
for i=1:N-1
for j=2:M-1
PH(i,j)=0.25*(phi(i+1,j+1)-phi(i+1,j-1)+phi(i,j+...
1)-phi(i,j-1));
UH(i,j)=g11H(i,j)*(phi(i+1,j)-phi(i,j))+g12H(i,j...
)*PH(i,j);
VH(i,j)=g21H(i,j)*(phi(i+1,j)-phi(i,j))+g22H(i,j...
)*PH(i,j);
end
end
% Calculamos la densidad, ecuación (4.13)
IMA=0;
for i=1:N
for j=1:M
DDV(i,j)=1-((dxdxiV(i,j)^2+dydxiV(i,j)^2)*UV(i,j)...
^2+(dxdetaV(i,j)^2+dydetaV(i,j)^2)*VV(i,j)^2+...
2*UV(i,j)*VV(i,j)*(dxdxiV(i,j)*dxdetaV(i,j)+...
dydxiV(i,j)*dydetaV(i,j)))/(2*H0);
DDH(i,j)=1-((dxdxiH(i,j)^2+dydxiH(i,j)^2)*UH(i,j)...
^2+(dxdetaH(i,j)^2+dydetaH(i,j)^2)*VH(i,j)^2+...
2*UH(i,j)*VH(i,j)*(dxdxiH(i,j)*dxdetaH(i,j)+...
dydxiH(i,j)*dydetaH(i,j)))/(2*H0);
if DDV(i,j)<0||DDH(i,j)<0
IMA=1;
end
dV(i,j)=d0*abs(DDV(i,j))^(1/(gamma-1));
dH(i,j)=d0*abs(DDH(i,j))^(1/(gamma-1));
end
end
% Introducimos las variables anteriores en la ecuación del potencial.
% Fórmula (4.11)
for i=1:N-1;
for j=2:M-1;
if i==1
phi(i,j)=(dH(i,j)*JH(i,j)*(g12H(i,j)*PH(i,j)+g11H(i,...
j)*phi(i+1,j))-dH(N-1,j)*JH(N-1,j)*(g12H(N-1,j)*...
PH(N-1,j)-g11H(N-1,j)*(phi(N-1,j)-C))+dV(i,j-1)...
*JV(i,j-1)*(g21V(i,j-1)*PV(i,j-1)+g22V(i,j-1)*...
phi(i,j))-dV(i,j)*JV(i,j)*(g21V(i,j)*PV(i,j)-...
g22V(i,j)*phi(i,j-1)))/(dH(i,j)*JH(i,j)*g11H(...
i,j)+dH(N-1,j)*JH(N-1,j)*g11H(N-1,j)+dV(i,j)*JV...
(i,j)*g22V(i,j)+dV(i,j-1)*JV(i,j-1)*g22V(i,j-1));
else
phi(i,j)=(dH(i,j)*JH(i,j)*(g12H(i,j)*PH(i,j)+g11H(i,...
j)*phi(i+1,j))-dH(i-1,j)*JH(i-1,j)*(g12H(i-1,j)*...
PH(i-1,j)-g11H(i-1,j)*(phi(i-1,j)))+dV(i,j-1)*JV(...
i,j-1)*(g21V(i,j-1)*PV(i,j-1)+g22V(i,j-1)*phi(i,...
j+1))-dV(i,j)*JV(i,j)*(g21V(i,j)*PV(i,j)-g22V(i,j...
)*phi(i,j-1)))/(dH(i,j)*JH(i,j)*g11H(i,j)+dH(i...
-1,j)*JH(i-1,j)*g11H(i-1,j)+dV(i,j)*JV(i,j)*g22V...
(i,j)+dV(i,j-1)*JV(i,j-1)*g22V(i,j-1));
end
% Aplicamos el método SOR de sobrerelajación, ecuación (4.29).
phi(i,j)=w*phi(i,j)+(1-w)*phiold(i,j);
end
end
%---------------------CONDICIÓN EN LA SUPEFICIE-----------------------%
% Aplicamos la fórmula (4.15)
for i=N-1:-1:2;
phi(i,M)=(1/3)*(4*phi(i,M-1)-phi(i,M-2)-g21(i,j)*(phi(i+...
1,M)-phi(i-1,M))/g22(i,j));
end
phi(1,M)=(1/3)*(4*phi(1,M-1)-phi(1,M-2)-g21(1,M)*(phi(2,M)...
-phi(N-1,M)+C)/g22(1,M));
%-------------CONDICIÓN EN LA DISCONTINUIDAD DEL POTENCIAL------------%
% Aplicamos la condición de Kutta, ecuación (4.16)
for j=1:M;
phi(N,j)=phi(1,j)+C;
end
ddd=max(max(abs(phi-phiold)));
%-----------------------CÁLCULO DE LA CIRCULACIÓN---------------------%
% Utilizamos la ecuación (4.18) que impone velocidad nula en el borde
% de salida
C=phi(N-1,M)-phi(2,M)-g12(1,M)*(phi(1,M-2)-4*phi(1,M-1)+...
3*phi(1,M))/g11(1,M);
end
end