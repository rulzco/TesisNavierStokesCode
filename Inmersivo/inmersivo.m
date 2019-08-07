%% --------------------------ESIME TICOMÁN---------------------------------
%%SOLVER Ehecatl-Versión 1.0.0 (Raúl Alberto Bernal Orozco, contacto:
%%rulozcohotmail.com. El presente solver es capaz de resolver las
%%ecuaciones de Navier-Stokes para un flujo incompresible y 2D sobre una
%%malla cartesiana,se han insertado funciones que permiten implementar un
%%método inmersivo de forzado directo, con lo que se simula el flujo
%%alrededor de un cilindro, sin embarga el código puede ser modificado para
%%resolver otros casos como cavidades, escalones, u otra geométrica que sea
%%debidamente implementada.
function [co,dt,nt,u,v,p,eu,ev,eP,bu,bv,iu,iv]=inmersivo
clc
tic
nx=300; % Número de nodos en el eje-x
ny=200; % Número de nodos en el eje-y
nit=10;
r=0.075; % Radio del cilindro [m]
L=2*r; 
xmin=0; xmax=3;
ymin=0; ymax=2;
dx=(xmax-xmin)/(nx-1);
dy=(ymax-ymin)/(ny-1);
x=linspace(0,xmax,nx);
y=linspace(0,ymax,ny);
[X,Y]=meshgrid(x,y);
t=0;% Tiempo inicial [s]
tmax=3.0;% Tiempo final [s]
dt=0.1*min(dx,dy);% Paso temporal [s] 
nt= round((tmax-t)/dt);
t0=ones(1,nt+1)*dt;

u_0=1; % Velocidad incial [m/s]
rho=1; % Densidad [kg/m^3]
nu=0.001; % Viscosidad cinemática [m^2/s]
Re=(u_0*L)/nu; % Número de Reynolds
co=u_0*dt/dx; % Número de Courant
%% --------------------Generación de la malla cartesiana------------------%
dis=linspace(0,pi,20);
xcirculo=r*cos(dis);
ycirculo=sqrt(r.^2-(xcirculo-0.0).^2);

xcirculo=[xcirculo,xcirculo(end:-1:1,:)];
ycirculo=[ycirculo,-ycirculo(end:-1:1,:)];
xcirculo=xcirculo+xmax*0.25;
ycirculo=ycirculo+ymax*0.5;
figure(1)
plot(X,Y,'g');hold on;
plot(X',Y','g');hold on;
plot(xcirculo,ycirculo,'k'); hold on;axis equal;
title('Malla Cartesiana');

iu=[]; % Vector que guarda los indices de los nodos dentro del cilindro
for i=1:length(x)
    for j=1:length(y)
        if sqrt((x(i)-xmax*0.25)^2+(y(j)-ymax*0.5)^2)>0 && sqrt((x(i)-xmax*0.25)^2+(y(j)-ymax*0.5)^2)<=r
            iu = cat(1, iu, [i,j]);
        end
    end
end


iv=[]; % Vector que guarda los indices de los nodos dentro del cilindro
for i=1:length(x)
    for j=1:length(y)
        if sqrt((x(i)-xmax*0.25)^2+(y(j)-ymax*0.5)^2)>0 && sqrt((x(i)-xmax*0.25)^2+(y(j)-ymax*0.5)^2)<=r
            iv = cat(1, iv, [i,j]);
        end
    end
end
            
bul=min(iu(:,1));
bur=max(iu(:,1));
bub=min(iu(:,2));
but=max(iu(:,2));

cux=linspace(bul,bur,bur-bul+1);

bu=[]; % Vector que guarda los indices de los nodos contiguos al cilindro

for i=1:length(cux)
    
    ind=find(iu(:,1)==cux(i));
 
    Min=iu(ind(1,1),2);
    Max=iu(ind(1,1),2);

    for j=1:length(ind)
        if iu(ind(j,1),2)<Min
            Min=iu(ind(j,1),2);
        end
        if iu(ind(j,1),2)>Max
            Max=iu(ind(j,1),2);
        end
    end
    bu = cat(1, bu, [cux(i),Min-1]);
    bu = cat(1, bu, [cux(i),Max+1]);
end

cuy=linspace(bub,but,but-bub+1);
for i=1:length(cuy)
    
    ind=find(iu(:,2)==cuy(i));
    
    Min=iu(ind(1,1),1);
    Max=iu(ind(1,1),1);
    
    for j=1:length(ind)
        if iu(ind(j,1),1)<Min
            Min=iu(ind(j,1),1);
        end
        if iu(ind(j,1),1)>Max
            Max=iu(ind(j,1),1);
        end
    end
    bu = cat(1, bu, [Min-1,cuy(i)]);
    bu = cat(1, bu, [Max+1,cuy(i)]);
end

figure(1)
for i=1:length(bu)
    plot(x(bu(i,1)),y(bu(i,2)),'.b'); hold on;
end

bvl=min(iv(:,1));
bvr=max(iv(:,1));
bvb=min(iv(:,2));
bvt=max(iv(:,2));

cvx=linspace(bvl,bvr,bvr-bvl+1);

bv=[]; % Vector que guarda los indices de los nodos contiguos al cilindro

for i=1:length(cvx)
    
    ind=find(iv(:,1)==cvx(i));
    
    Min=iv(ind(1,1),2);
    Max=iv(ind(1,1),2);

    for j=1:length(ind)
        if iv(ind(j,1),2)<Min
            Min=iv(ind(j,1),2);
        end
        if iv(ind(j,1),2)>Max
            Max=iv(ind(j,1),2);
        end
    end
    % bv=[cvx(i),Min-1;cvx(i),Max+1];
    bv = cat(1, bv, [cvx(i),Min-1]);
    bv = cat(1, bv, [cvx(i),Max+1]);
end

cvy=linspace(bvb,bvt,bvt-bvb+1);

for i=1:length(cvy)
    
    ind=find(iv(:,2)==cvy(i));
    
    Min=iv(ind(1,1),1);
    Max=iv(ind(1,1),1);
    
    for j=1:length(ind)
        if iv(ind(j,1),1)<Min
            Min=iv(ind(j,1),1);
        end
        if iv(ind(j,1),1)>Max
            Max=iv(ind(j,1),1);
        end
    end
    % bv=[bv;Min-1,cvy(i);Max+1,cvy(i)];
    bv = cat(1, bv, [Min-1,cvy(i)]);
    bv = cat(1, bv, [Max+1,cvy(i)]);
end

figure(1)
for i=1:length(bv)
    plot(x(bv(i,1)),y(bv(i,2)),'xk');hold on;
    xlabel('L [m]');ylabel('h [m]')
end

%% Inicializa las matrices
us=zeros(nx,ny); ud=zeros(nx,ny); udd=zeros(nx,ny);
vs=zeros(nx,ny); vd=zeros(nx,ny); vdd=zeros(nx,ny);
ps=zeros(nx,ny);
%b=zeros(nx,ny);

%% -------------------------Condiciones de frontera-----------------------%
[us,vs]=condicionesFrontera(us,vs);

p_prime = zeros(nx,ny); A_star= zeros(nx-1,ny-1);
p_prime2= zeros(nx,ny); B_star= zeros(nx-1,ny-1);
u_prime = zeros(nx,ny); d= zeros(nx-1,ny-1);
v_prime = zeros(nx,ny);

p = ps+p_prime;
u = us+u_prime;
v = vs+v_prime;

%% --------------------Solver de las ecuaciones de Navier-Stokes----------%
% -------------------Ecuaciones () a ()-------------------------------%
alpha = 0.4; 
h = waitbar(0,'Resolviendo N-S...');k=1;
for itr = 1:nt+1
    un=u;
    vn=v;
    for i = 2:nx-1
        for j = 2:ny-1
            udd(i,j) = 0.5 * (u(i-1,j) + u(i-1,j+1));
            ud(i,j)  = 0.5 * (u(i,j) + u(i,j+1));
            vd(i,j)  = 0.5 * (v(i,j) + v(i+1,j));
            vdd(i,j) = 0.5 * (v(i,j-1) + v(i+1,j-1));
            A_star(i,j) =  -((u(i,j+1)*vd(i,j) - u(i,j-1)*vdd(i,j))/(2*dy) + ((u(i+1,j))^2 - (u(i-1,j))^2)/(2*dx))...
                +nu*((u(i+1,j)-2*u(i,j)+u(i-1,j))/dx^2 + (u(i,j+1)-2*u(i,j)+u(i,j-1))/dy^2);
            B_star(i,j) =  -((v(i+1,j)*ud(i,j) - v(i-1,j)*udd(i,j))/(2*dx) + ((v(i,j+1))^2 - (v(i,j-1))^2)/(2*dy))...
                +nu*((v(i+1,j)-2*v(i,j)+v(i-1,j))/dx^2 + (v(i,j+1)-2*v(i,j)+v(i,j-1))/dy^2);
            us(i,j) = us(i,j) + A_star(i,j)*dt - dt/dx*(ps(i+1,j) - ps(i,j))/rho;
            vs(i,j) = vs(i,j) + B_star(i,j)*dt - dt/dy*(ps(i,j+1) - ps(i,j))/rho;
        end
    end
    [us,vs]=condicionesFrontera(us,vs); %Condiciones de fontera velocidad

    a=2*( dt/dx^2 + dt/dy^2 );
	b=-dt/dx^2;
	c=-dt/dy^2;
    lambda=0.1;
    
    for itrate = 1:nit
        for i = 2:nx-1
            for j = 2:ny-1
                % p_prime2(i,j) = 0.25*(p_prime(i-1,j)+p_prime(i+1,j)+p_prime(i,j-1)+p_prime(i,j+1)) + dt*(us(i,j)-us(i+1,j)+vs(i,j)-vs(i,j+1));
                d(i,j)=(1/dx)*(us(i,j) - us(i-1,j)) + (1/dy)*( vs(i,j) - vs(i,j-1));
                p_prime2(i,j) = -(1/a)*( b*p_prime(i+1,j) + b*p_prime(i-1,j) + c*p_prime(i,j+1) + c*p_prime(i,j-1) + d(i,j));
                p_prime2(i,j)=lambda*p_prime2(i,j)+(1-lambda)*p_prime(i,j);  
            end
        end
        p_prime2(:,ny)=p_prime2(:,ny-1); % dp/dy=o @ y=0
        p_prime2(1,:)=0;     % dp/dx=o @ x=0
        p_prime2(:,1)=p_prime2(:,2);     % dp/dx=o @ x=0
        p_prime2(nx,:)=0;
        % p_prime2(:,ny)=0;         % dp=0 @ y=2
        mx = max(max(abs(p_prime-p_prime2)));
        p_prime = p_prime2;
        if mx<0.001, break, end
    end
    
    for j=2:(ny-1)
        for i=2:(nx-1)
            v(i,j) = vs(i,j) - (dt/dy) * (p_prime(i,j+1) - p_prime(i,j));
			vs(i,j)=v(i,j);
	
			u(i,j) = us(i,j) - (dt/dx) * (p_prime(i+1,j) - p_prime(i,j));
			us(i,j)=u(i,j);
            
        end
    end
    
    p = ps + alpha*p_prime;
    ps = p;
    p_prime=zeros(nx,ny);
    
% ------Aplica los terminos del forzado directo en el interior y la
% frontera del cilindro---------------------------------------------------%
    for i=1:length(bu)
        us(bu(i,1),bu(i,2))=0;
        %vs(bv(i,1),bv(i,2))=0;
        % p(bu(i,1),bu(i,2))=p(bu(i,1)+1,bu(i,2)+1);
    end
    for i=1:length(iu)
        us(iu(i,1),iu(i,2))=0;
        p(iu(i,1),iu(i,2))=0;
    end
    for i=1:length(bv)
        vs(bv(i,2),bv(i,1))=0;
        % p(bv(i,1),bv(i,2))=p(bu(i,1)+1,bv(i,2)+1);
    end
    for i=1:length(iv)
        vs(iv(i,1),iv(i,2))=0;
        p(iv(i,1),iv(i,2))=0;
    end
% ------------------------------------------------------------------------%    
    ev(1,itr)=max(max(abs(v-vn)))/max(max(abs(v)));
    eu(1,itr)=max(max(abs(u-un)))/max(max(abs(u)));
    eP(1,itr)=mx; 
    %t = cumsum(t0);
    t(1,itr) = 0+itr*dt;
    
    waitbar(itr / (nt+1))
    
    if k==1;
        plotFields(x,y,p,u,v,xcirculo,ycirculo,xmax,ymax,itr,dt);
    elseif k==25; 
        plotFields(x,y,p,u,v,xcirculo,ycirculo,xmax,ymax,itr,dt);
        k=0;
    end
    k=k+1;
    
    for i=1:length(iu)
        u(iu(i,1),iu(i,2))=0;
    end

    for i=1:length(iv)
        v(iv(i,1),iv(i,2))=0;
    end
    
    save u.mat u;
    save v.mat v;
    save p.mat p;
    save eu.mat eu;
    save ev.mat ev;
    save ep.mat eP;
    save t.mat t;
    
    if itr==1
       q=figure
        set(q,'WindowStyle','docked');
    end
    figure(6); clf
    q1=semilogy(t,ev,'r'); hold on;
    q2=semilogy(t,eu,'b'); hold on;
    q3=semilogy(t,eP,'g'); hold on; xlabel('Tiempo [s]');ylabel('Error');
    title('Residuales');legend([q1 q2 q3],'Velocidad v','Velocidad u','Presión');
      
end
delete(h) 
toc
%{
figure(2)
surf(x,y,p','Facealpha',1,'EdgeColor','none');view(0,90);hold on; % Campo presión
colormap(jet); % parula jet hsv cool hot spring summer autmn winter
shading interp
% quiver(x,y,u',v',2,'k');hold on; % grafica vectores de velocidad
fill3(xcirculo,ycirculo,ycirculo*0+abs(max(max(p)))*0.5e3,'k');hold on;
axis image,colorbar;
xlabel('x'); ylabel('y');
title({'Campo de Presión';['{Tiempo} = ',num2str( itr(end)*dt ),'s']});

figure(3)
surf(x,y,sqrt(u'.^(2)+v'.^(2)),'Facealpha',1,'EdgeColor','flat');view(0,90);hold on;% Campo de velocidad
%contourf(x,y,sqrt(u'.^(2)+v'.^(2)),'linewidth',0);
colormap(jet);% parula jet hsv cool hot spring summer autmn winter
shading interp
axis image,colorbar;
%patch(xcirculo,ycirculo,'k');hold on;
fill3(xcirculo,ycirculo,ycirculo*0+max(max(u.^2+v.^2)*5),'k');hold on;
title({'Campo de velocidad [m/s]';['{Tiempo} = ',num2str( itr(end)*dt ),'s']});

figure(5);clf;
streamslice(x,y,u',v',2);hold on;
fill(xcirculo,ycirculo,'k');
axis equal;
xlabel('X');ylabel('Y')
title({'Líneas de corriente';['{Tiempo} = ',num2str( itr(end)*dt ),'s']});
drawnow
frame5 = getframe(5);
im5{itr} = frame2im(frame5);

figure(6); clf
q1=semilogy(t,ev,'r'); hold on;
q2=semilogy(t,eu,'b'); hold on;
q3=semilogy(t,eP,'g'); hold on; xlabel('Tiempo [s]');ylabel('Error');
title('Residuales');legend([q1 q2 q3],'Velocidad v','Velocidad u','Presión');
%}
%{
load u.mat;
load v.mat;
load p.mat;
load eu.mat;
load ev.mat;
load ep.mat;
load t.mat;
%}


fprintf(' Diámetro del cilindro %.4f \n Re %.4f \n Courant %.5f \n dt %.5f ',L,Re,co,dt);
fprintf(' Error U %.5f \n Error V %.5f \n Error P %.5f ',eu(length(eu))...
    ,ev(length(ev)),eP(length(eP)));

    function [us,vs]=condicionesFrontera(us,vs)
        us(1,:)=u_0;
        us(:,1)=u_0; % Establece la velocidad igual a la velocidad incial
        us(nx,:)=u_0;%us(nx-1,:);
        us(:,ny)=u_0;   
        vs(1,:)=0;
        vs(:,ny)=0;
        vs(:,1)=0;
        vs(nx,:)=0;%vs(nx-1,:);
    end

    function plotFields(x,y,p,u,v,xcirculo,ycirculo,xmax,ymax,itr,dt)
        figure(2)
        surf(x,y,p','Facealpha',1,'EdgeColor','none');view(0,90);hold on; % Campo presión
        colormap(jet); % parula jet hsv cool hot spring summer autmn winter
        shading interp;%caxis([min(min(p)) max(max(p))]);
        % quiver(x,y,u',v',2,'k');hold on; % grafica vectores de velocidad
        fill3(xcirculo,ycirculo,ycirculo*0+abs(max(max(p)))*0.5e3,'k');hold on;
        axis image,colorbar;
        xlabel('x'); ylabel('y');
        title({'Campo de Presión';['{Tiempo} = ',num2str( itr*dt ),'s']});
        drawnow
        frame2 = getframe(2);
        im2{itr} = frame2im(frame2);
        
        filename = 'Campo de Presión.gif'; % Específica el tipo de archivo generado
        [A,map] = rgb2ind(im2{itr},256);
        if itr == 1
            imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.3);
        else
            imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.3);
        end
        
        figure(3)
        surf(x,y,sqrt(u'.^(2)+v'.^(2)),'Facealpha',1,'EdgeColor','flat');view(0,90);hold on;% Campo de velocidad
        %contourf(x,y,sqrt(u'.^(2)+v'.^(2)),'linewidth',0);
        colormap(jet);% parula jet hsv cool hot spring summer autmn winter
        shading interp
        axis image,colorbar;
        %patch(xcirculo,ycirculo,'k');hold on;
        fill3(xcirculo,ycirculo,ycirculo*0+max(max(u.^2+v.^2)*5),'k');hold on;
        title({'Campo de velocidad [m/s]';['{Tiempo} = ',num2str( itr*dt ),'s']});
        drawnow
        frame3 = getframe(3);
        im3{itr} = frame2im(frame3);
        
        filename = 'Campo de velocidad.gif'; % Específica el tipo de archivo generado
        [A,map] = rgb2ind(im3{itr},256);
        if itr == 1
            imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.3);
        else
            imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.3);
        end
        
        figure(4);clf;
        quiver(x,y,u',v',2,'k');hold on;% Gráfica vectores de velocidad
        fill(xcirculo,ycirculo,'k');
        axis equal; 
        xlabel('x'); ylabel('y');
        title({'Vectores de Velocidad';['{Tiempo} = ',num2str( itr*dt ),'s']});
        drawnow
        frame4 = getframe(4);
        im4{itr} = frame2im(frame4);
        
        filename = 'Vectores de Velocidad.gif'; % Específica el tipo de archivo generado
        [A,map] = rgb2ind(im4{itr},256);
        if itr == 1
            imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.3);
        else
            imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.3);
        end
        
        figure(5);clf;
        streamslice(x,y,u',v',2);hold on;
        fill(xcirculo,ycirculo,'k');
        axis equal;
        xlabel('X');ylabel('Y')
        title({'Líneas de corriente';['{Tiempo} = ',num2str( itr*dt ),'s']});
        drawnow
        frame5 = getframe(5);
        im5{itr} = frame2im(frame5);
        
        filename = 'Líneas de corriente.gif'; % Específica el tipo de archivo generado
        [A,map] = rgb2ind(im5{itr},256);
        if itr == 1
            imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.3);
        else
            imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.3);
        end
    end
end