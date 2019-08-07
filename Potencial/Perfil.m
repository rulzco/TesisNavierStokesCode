function [x,y] = Perfil(cN,m,p,t,N)
%----------------------------VARIABLES--------------------------------%
%thetan= ángulo que forman las tangentes a la línea de curvatura media
%Xn= distribución de la variable x sobre la cuerda
%Yc= posición de la línea media
%Yt= espesor con respecto a la línea media
%x,y= vectores de posición de la superficie
%---------------------------------------------------------------------%
% Inicialmente vamos a considerar un valor de la cuerda unidad.
m=m/100;p=p/10;t=t/100;
% Calculamos la dimensión del vector horizontal Xn. Como en el borde
% de salida tenemos dos líneas y en el borde de ataque sólo una, N
% tiene que ser un número impar.
Nnaca=(N+1)/2;
thetaXn=linspace(0,1,Nnaca);
Xn=0.5*(1-cos(pi.*thetaXn));
Yc=zeros(Nnaca,1);
Yt=zeros(Nnaca,1);
thetan=zeros(Nnaca,1);
% Aplicamos las fórmulas dadas en (1.1) hasta (1.10).
    for i=1:Nnaca
        if Xn(i)<p
            Yc(i)=(m/(p^2))*(2*p*Xn(i)-Xn(i)^2);
            thetan(i)=atan(((2*m)/p^2)*(p-Xn(i)));
        else
            Yc(i)=(m/(1-p)^2)*(1-2*p+2*p*Xn(i)-Xn(i)^2);
            thetan(i)=atan((((2*m)/(1-p)^2))*(p-Xn(i)));
        end
        Yt(i)=(t/0.2)*(0.2969*sqrt(Xn(i))-0.126*Xn(i)- 0.3516*Xn(i)^2+...
        0.2843*Xn(i)^3-0.1015*Xn(i)^4);
    end
    YTs=Yt.*sin(thetan);
    YTc=Yt.*cos(thetan);
% Almacenamos los valores de posición según nos encontremos en el
% extradós o intradós del perfil.
    for j=Nnaca-1:-1:1
        i=Nnaca+1-j;
        x(i)=Xn(j)-YTs(j);
        y(i)=Yc(j)+YTc(j);
    end
        x(1)=1;
        y(1)=0;
        x(N)=x(1);
        y(N)=y(1);
    for j=2:Nnaca-1
        i=Nnaca+j-1;
        x(i)=Xn(j)+YTs(j);
        y(i)=Yc(j)-YTc(j);
    end
% Escalamos los valores x e y según el valor de la cuerda.
x=x*cN;y=y*cN;
end