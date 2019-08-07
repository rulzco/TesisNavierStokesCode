function [X,Y] = Malla(N,M,x,y)
%-----------------------------VARIABLES-------------------------------%
%xi,eta= Coordenadas arbitrarias
%X,Y= Coordenadas en los planos x,y
%ite= Número de iteracción
%itemax= Número de iteracción máxima que imponemos
%ddd1= Diferencia de la coordenada x con respecto a la iteracción anterior
%ddd2= Diferencia de la coordenada y con respecto a la iteracción anterior
%tol= Tolerancia admitida
%Q= Función para generar atracción en una línea coordenada
%etak, Ak, Ck= Parámetros que regulan la atracción de la línea coordenada
%w= Parámetro que regula la velocidad de convergencia de la solución
%---------------------------------------------------------------------%
%-------------------CONDICIÓN SUPERFICIE DEL PERFIL-------------------%
% Almacenamos en la superficie el valor obtenido de la función perfil.
X(:,M)=x(:);
Y(:,M)=y(:);
%---------------------CONDICIÓN LEJOS DEL PERFIL----------------------%
R=20;
theta=linspace(0,2*pi,N);
X(:,1)=0.5+R*cos(theta);
Y(:,1)=R*sin(theta);
X=X-0.5;
%--------------------CONDICIÓN INICIAL ARBITRARIA---------------------%
for j=2:M-1
X(:,j)=X(:,j-1)+(X(:,M)-X(:,1))/(M-1);
Y(:,j)=Y(:,j-1)+(Y(:,M)-Y(:,1))/(M-1);
end
%-----------------DESARROLLO DE LA FUNCIÓN DE LAPLACE-----------------%
etak=M;Ak=20;Ck=0.2;
eta=1:M;xi=1:N;
ite=0;itemax=8000;tol=1.e-7;ddd1=1;ddd2=1;w=1.8;
while ddd1>tol&&ddd2>tol&&ite<itemax
ite=ite+1;
Xold=X;
Yold=Y;
% Aumento de la densidad del mallado en la superficie
% Fórmula (1.2)
for i=1:N
for j=1:M
Q(i,j)=Ak*exp(-Ck*abs(eta(j)-etak));
P(i,j)=0;
end
end
%---------------CONDICIÓN NODOS INTERIORES DE LA MALLA----------------%
% Desarrollo de las ecuaciones (4.1) hasta la (4.6)
for j=2:M-1
for i=2:N
if i==N
alfa(N,j)=0.25*((X(i,j+1)-X(i,j-1))^2+(Y(i,j+1)-Y(...
i,j-1))^2);
beta(N,j)=0.25*((X(i,j+1)-X(i,j-1))*(X(2,j)-X(i-1,...
j)) +(Y(2,j)-Y(i-1,j))*(Y(i,j+1)-Y(i,j-1)));
gammaM(N,j)=0.25*((X(2,j)-X(i-1,j))^2+(Y(2,j)-Y(i-...
1,j))^2);
J(N,j)=0.25*((X(2,j)-X(i-1,j))*(Y(i,j+1)-Y(i,j-1))...
-(Y(2,j)-Y(i-1,j))*(X(i,j+1)-X(i,j-1)));
X(i,j)=(0.5/(alfa(i,j)+gammaM(i,j)))*(alfa(i,j)*(X(...
2,j)+X(i-1,j))+gammaM(i,j)*(X(i,j+1)+X(i,j-1))-...
0.5*beta(i,j)*(X(2,j+1)-X(2,j-1)-X(i-1,j+1)+X(i...
-1,j-1))+Q(i,j)*0.5*J(i,j)^2*(X(i,j+1)-X(i,j-1)...
)+P(i,j)*0.5*J(i,j)^2*(X(2,j)-X(i-1,j)));
Y(i,j)=(0.5/(alfa(i,j)+gammaM(i,j)))*(alfa(i,j)*(Y(...
2,j)+Y(i-1,j))+gammaM(i,j)*(Y(i,j+1)+Y(i,j-1))-...
0.5*beta(i,j)*(Y(2,j+1)-Y(2,j-1)-Y(i-1,j+1)+Y(i...
-1,j-1))+Q(i,j)*0.5*J(i,j)^2*(Y(i,j+1)-Y(i,j-1)...
)+P(i,j)*0.5*J(i,j)^2*(Y(2,j)-Y(i-1,j)));
else
alfa(i,j)=0.25*((X(i,j+1)-X(i,j-1))^2+(Y(i,j+1)-Y(i...
,j-1))^2);
beta(i,j)=0.25*((X(i,j+1)-X(i,j-1))*(X(i+1,j)-X(i-...
1,j))+(Y(i+1,j)-Y(i-1,j))*(Y(i,j+1)-Y(i,j-1)));
gammaM(i,j)=0.25*((X(i+1,j)-X(i-1,j))^2+(Y(i+1,j)-...
Y(i-1,j))^2);
J(i,j)=0.25*((X(i+1,j)-X(i-1,j))*(Y(i,j+1)-Y(i,j-1)...
)-(Y(i+1,j)-Y(i-1,j))*(X(i,j+1)-X(i,j-1)));
X(i,j)=(0.5/(alfa(i,j)+gammaM(i,j)))*(alfa(i,j)*(X(...
i+1,j)+X(i-1,j))+gammaM(i,j)*(X(i,j+1)+X(i,j-1)...
)-0.5*beta(i,j)*(X(i+1,j+1)-X(i+1,j-1)-X(i-1,j+...
1)+X(i-1,j-1))+Q(i,j)*0.5*J(i,j)^2*(X(i,j+1)-X(...
i,j-1))+P(i,j)*0.5*J(i,j)^2*(X(i+1,j)-X(i-1,j)));
Y(i,j)=(0.5/(alfa(i,j)+gammaM(i,j)))*(alfa(i,j)*(Y(...
i+1,j)+Y(i-1,j))+gammaM(i,j)*(Y(i,j+1)+Y(i,j-1)...
)-0.5*beta(i,j)*(Y(i+1,j+1)-Y(i+1,j-1)-Y(i-1,j+...
1)+Y(i-1,j-1))+Q(i,j)*0.5*J(i,j)^2*(Y(i,j+1)-Y(...
i,j-1))+P(i,j)*0.5*J(i,j)^2*(Y(i+1,j)-Y(i-1,j)));
end
% Aplicamos el método SOR de sobrerelajación, ecuación (4.29).
X(i,j)=w*X(i,j)+(1-w)*Xold(i,j);
Y(i,j)=w*Y(i,j)+(1-w)*Yold(i,j);
end
end
X(1,:)=X(N,:);
Y(1,:)=Y(N,:);
ddd1=max(max(abs(X-Xold)));
ddd2=max(max(abs(Y-Yold)));
end
X=X+0.5;
end