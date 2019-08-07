function [u,v] = Velocidades(alf,C,Machinf,theta,X,Y,phi,N,M,Vinf)
[g11,g22,g21,g12,J,dxdxi,dxdeta,dydxi,dydeta]=Tensor(X,Y,N,M);
%--------------------CONDICIÓN FRONTERA EXTERIOR----------------------%
% Ecuaciones (4.21) y (4.22).
    alfa=alf*pi/180;
    for i=1:N
        j=1;
        u(i,j)=Vinf*cos(alfa)+(C/(2*pi))*sqrt(1-Machinf^2)*(1+tan(...
        theta(i,j)-alfa)^2)*(1/(1+(Y(i,j)/X(i,j))^2))*(-Y(i,j)/...
        (X(i,j))^2)/(1+(1-Machinf^2)*(tan(theta(i,j)-alfa))^2);
        v(i,j)=Vinf*sin(alfa)+(C/(2*pi))*sqrt(1-Machinf^2)*(1+tan(...
        theta(i,j)-alfa)^2)*(1/(1+(Y(i,j)/X(i,j))^2))*(1/X(i,j)...
        )/(1+(1-Machinf^2)*(tan(theta(i,j)-alfa))^2);
    end
    %--------------------CONDICIÓN NODOS INTERIORES-----------------------%
    % Ecuaciones (4.23) y (4.24).
    for i=2:N-1
        for j=2:M-1
             u(i,j)=(1/J(i,j))*(((phi(i+1,j)-phi(i-1,j))/2)*dydeta(i,j...
             )-((phi(i,j+1)-phi(i,j-1))/2)*dydxi(i,j));
             v(i,j)=(1/J(i,j))*(((phi(i,j+1)-phi(i,j-1))/2)*dxdxi(i,j)...
             -((phi(i+1,j)-phi(i-1,j))/2)*dxdeta(i,j));
        end
    end
    for j=2:M-1
        u(1,j)=(1/J(1,j))*(((phi(2,j)-phi(N-1,j)+C)/2)*dydeta(1,...
        j)-((phi(1,j+1)-phi(1,j-1))/2)*dydxi(1,j));
        v(1,j)=(1/J(1,j))*(((phi(1,j+1)-phi(1,j-1))/2)*dxdxi(1,j)-((...
        phi(2,j)-phi(N-1,j)+C)/2)*dxdeta(1,j));
    end
    u(N,:)=u(1,:);
    v(N,:)=v(1,:);
    %--------------------CONDICIÓN FRONTERA INTERIOR----------------------%
    % Ecuaciones (4.26) y (4.27).
    j=M;
    for i=2:N-1
        u(i,j)=(((phi(i+1,j)-phi(i-1,j))/2)/J(i,j))*(dydeta(i,j)+...
        dydxi(i,j)*g21(i,j)/g22(i,j));
        v(i,j)=-(((phi(i+1,j)-phi(i-1,j))/2)/J(i,j))*(dxdeta(i,j)+...
        dxdxi(i,j)*g21(i,j)/g22(i,j));
    end
end