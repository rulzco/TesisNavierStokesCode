function [psi,Mach] = Corriente(u,v,gamma,H0,d0,p,X,Y,N,M)
% Inicializamos a cero la matriz corriente.
psi=zeros(N,M);
[g11,g22,g21,g12,J,dxdxi,dxdeta,dydxi,dydeta]=Tensor(X,Y,N,M);
d=d0.*(1-(u.^2+v.^2)/(2*H0)).^(1/(gamma-1));
    for i=1:N
        for j=1:M
            JU(i,j)=u(i,j)*dydeta(i,j)-v(i,j)*dxdeta(i,j);
        end
    end
% Como la función corriente es nula en la superficie, no modificamos
% su valor.
% Aplicamos la ecuación (4.19).
    for i=1:N
        for j=M:-1:2
            psi(i,j-1)=psi(i,j)+JU(i,j)*d(i,j);
        end
    end
    for i=1:N
        for j=1:M
            c(i,j)=sqrt(gamma*p(i,j)/d(i,j));
            Mach(i,j)=sqrt(u(i,j)^2+v(i,j)^2)/c(i,j);
        end
    end
end