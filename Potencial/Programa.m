clc; clear;
%-----------------------------VARIABLES-------------------------------%
%phi= Funci�n potencial
%psi= Funci�n de corriente
%u,v= Vectores velocidad seg�n las direcciones x e y
%Cp= Coeficiente de presi�n
%C= Circulaci�n alrededor del perfil
%Mach= N�mero de Mach
%L= Sustentaci�n
%D= Resistencia
%N,M= N�mero de l�neas coordenadas en xi,eta respectivamente
%alf= �ngulo de ataque
%cN= Cuerda del perfil
%m, p, t= par�metros caracter�sticos del perfil NACA
%Vinf= M�dulo de la velocidad lejos del perfil
%Tinf= Temperatura lejos del perfil
%pinf= Presi�n lejos del perfil
%gamma= Coeficiente adiab�tico del aire
%cp= Capacidad calor�fica a presi�n constante
%c= Velocidad del sonido
% Los sub�ndices 0 significan condici�n de remanso
% Los sub�ndices inf significan condici�n lejos del perfil
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
%--------------------------FUNCI�N PERFIL-----------------------------%
[x,y] = Perfil(cN,m,p,t,N);
%--------------------Representaci�n del perfil------------------------%
figure(1); clf
plot(x(:),y(:),'k');
%--------------------------FUNCI�N MALLA------------------------------%
[X,Y] = Malla(N,M,x,y);
%--------------------Representaci�n de la malla-----------------------%
figure(2); clf
for j=1:M
plot(X(:,j),Y(:,j)); hold on;
end
for i=1:N
plot(X(i,:),Y(i,:)); hold on;
end
[g11,g22,g21,g12,J,dxdxi,dxdeta,dydxi,dydeta]=Tensor(X,Y,N,M);
%------------------------FUNCI�N POTENCIAL----------------------------%
[phi,C,theta,IMA] = Potencial(d0,H0,gamma,Machinf,X,Y,Vinf,alf,N,M);
%----------------Representaci�n l�neas equipotenciales----------------%
figure(3); clf
plot(X(:,M),Y(:,M),'k');
hold on
xlabel('X');ylabel('Y');
contour(X,Y,phi,100);
%-----------------------FUNCI�N VELOCIDADES---------------------------%
[u,v] = Velocidades(alf,C,Machinf,theta,X,Y,phi,N,M,Vinf);
%----------------Representaci�n vectores de velocidad-----------------%
figure(4); clf
plot(X(:,M),Y(:,M),'k');
hold on
xlabel('X');ylabel('Y');
quiver(X,Y,u,v);
colorbar();
%-----------------FUNCI�N DISTRIBUCI�N DE PRESION---------------------%
[Cp,p]=Presiones(u,v,Vinf,dinf,gamma,pinf,p0,d0,H0);
%-------------------Representaci�n l�neas is�baras--------------------%
figure(5); clf
plot(X(:,M),Y(:,M),'k');
hold on
xlabel('X');ylabel('Y');
contour(X,Y,Cp,30);
%--------------Representaci�n presi�n sobre la superficie-------------%
figure(6); clf
xlabel('X');ylabel('Cp');
plot(X(1:(N-1)/2,M),Cp(1:(N-1)/2,M),'r');
hold on
plot(X((N-1)/2:N,M),Cp((N-1)/2:N,M),'b');
%-----------------------FUNCI�N DE CORRIENTE--------------------------%
[psi,Mach] = Corriente(u,v,gamma,H0,d0,p,X,Y,N,M);
%-----------------Representaci�n l�neas de corriente------------------%
figure(7); clf
plot(X(:,M),Y(:,M),'k');
hold on
xlabel('X');ylabel('Y');
psin=real(psi);Machn=real(Mach);
contour(X,Y,psin,50);
%------------------Representaci�n contornos de Mach-------------------%
figure(8); clf
plot(X(:,M),Y(:,M),'k');
hold on
xlabel('X');ylabel('Y');
contour(X,Y,Machn,50);
%----------------FUNCI�N COEFICIENTES AERODIN�MICOS-------------------%
[L,D]=Coeficientes(N,M,Cp,X,Y,alf,cN);
if IMA==1
disp( 'No converge')
end
