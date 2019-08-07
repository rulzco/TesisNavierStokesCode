clc; clear;
%-----------------------------VARIABLES-------------------------------%
%phi= Función potencial
%psi= Función de corriente
%u,v= Vectores velocidad según las direcciones x e y
%Cp= Coeficiente de presión
%C= Circulación alrededor del perfil
%Mach= Número de Mach
%L= Sustentación
%D= Resistencia
%N,M= Número de líneas coordenadas en xi,eta respectivamente
%alf= Ángulo de ataque
%cN= Cuerda del perfil
%m, p, t= parámetros característicos del perfil NACA
%Vinf= Módulo de la velocidad lejos del perfil
%Tinf= Temperatura lejos del perfil
%pinf= Presión lejos del perfil
%gamma= Coeficiente adiabático del aire
%cp= Capacidad calorífica a presión constante
%c= Velocidad del sonido
% Los subíndices 0 significan condición de remanso
% Los subíndices inf significan condición lejos del perfil
%---------------------------------------------------------------------%
Tinf=273.15;
pinf=101325;
Vinf=75;
d=1;
cN=d;
N=35;%division sobre el eje xi
M=35;%division sobre el eje eta
m=1;p=2;t=10;cN=1;
alf=1;
% Calculamos las condiciones del aire
gamma=1.4;
cp=1006;
Rg=cp*(gamma-1)/gamma;
dinf=pinf/(Rg*Tinf);
Hinf=cp*Tinf;
cinf=sqrt(gamma*pinf/dinf);
% Calculamos las condiciones de remanso del aire
H0=Hinf+0.5*Vinf^2;
d0=dinf/(1-0.5*Vinf^2/H0);
p0=pinf*(d0/dinf)^gamma;
% Comprobamos inicialmente que las condiciones del fluido se encuentren
% dentro de las condiciones impuestas.
Machinf=Vinf/cinf;
Re=Vinf*cN*dinf/17e-6;

if Machinf>0.8||Re>10^7
disp('El flujo no cumple las condiciones necesarias')
exit
end
%--------------------------FUNCIÓN PERFIL-----------------------------%
[x,y] = Perfil(cN,m,p,t,N);
%--------------------Representación del perfil------------------------%
figure(1); clf
plot(x(:),y(:),'k');
%--------------------------FUNCIÓN MALLA------------------------------%
[X,Y] = Malla(N,M,x,y);
%--------------------Representación de la malla-----------------------%
figure(2); clf
for j=1:M
plot(X(:,j),Y(:,j)); hold on;
end
for i=1:N
plot(X(i,:),Y(i,:)); hold on;
end
[g11,g22,g21,g12,J,dxdxi,dxdeta,dydxi,dydeta]=Tensor(X,Y,N,M);
%------------------------FUNCIÓN POTENCIAL----------------------------%
[phi,C,theta,IMA] = Potencial(d0,H0,gamma,Machinf,X,Y,Vinf,alf,N,M);
%----------------Representación líneas equipotenciales----------------%
figure(3); clf
plot(X(:,M),Y(:,M),'k');
hold on
xlabel('X');ylabel('Y');
contour(X,Y,phi,100);
%-----------------------FUNCIÓN VELOCIDADES---------------------------%
[u,v] = Velocidades(alf,C,Machinf,theta,X,Y,phi,N,M,Vinf);
%----------------Representación vectores de velocidad-----------------%
figure(4); clf
plot(X(:,M),Y(:,M),'k');
hold on
xlabel('X');ylabel('Y');
quiver(X,Y,u,v);
colorbar();
%-----------------FUNCIÓN DISTRIBUCIÓN DE PRESION---------------------%
[Cp,p]=Presiones(u,v,Vinf,dinf,gamma,pinf,p0,d0,H0);
%-------------------Representación líneas isóbaras--------------------%
figure(5); clf
plot(X(:,M),Y(:,M),'k');
hold on
xlabel('X');ylabel('Y');
contour(X,Y,Cp,30);
%--------------Representación presión sobre la superficie-------------%
figure(6); clf
xlabel('X');ylabel('Cp');
plot(X(1:(N-1)/2,M),Cp(1:(N-1)/2,M),'r');
hold on
plot(X((N-1)/2:N,M),Cp((N-1)/2:N,M),'b');
%-----------------------FUNCIÓN DE CORRIENTE--------------------------%
[psi,Mach] = Corriente(u,v,gamma,H0,d0,p,X,Y,N,M);
%-----------------Representación líneas de corriente------------------%
figure(7); clf
plot(X(:,M),Y(:,M),'k');
hold on
xlabel('X');ylabel('Y');
psin=real(psi);Machn=real(Mach);
contour(X,Y,psin,50);
%------------------Representación contornos de Mach-------------------%
figure(8); clf
plot(X(:,M),Y(:,M),'k');
hold on
xlabel('X');ylabel('Y');
contour(X,Y,Machn,50);
%----------------FUNCIÓN COEFICIENTES AERODINÁMICOS-------------------%
[L,D]=Coeficientes(N,M,Cp,X,Y,alf,cN);
if IMA==1
disp( 'No converge')
end
