function [t,s,M]=waveation()
% 2D Shallow Water Model
% Lax-Wendroff finite difference method.
% [t,s,M]=waveation()
% t = time intervals
% s = energy of the wave at certain distance for each time interval
% M = frames of the simulation

% Parameters

n=100;                  % grid size
g=9.8;                 % gravitational constant
dt=0.001;               % hardwired timestep
dx=0.1;                 % space step
dy=0.1;
nplotstep=10;           % plot interval
nn=0.5;                 % viscous friction factor

x=1:n+2;
y=1:n+2;
[yy,xx]=meshgrid(y,x);

%topography

%almost flat

Z=yy/(10*(n+2));

% steep coast

% Z=0.9*(yy-(n-8))/10;
% Z(Z<0)=0;

%boundary conditions

%simple

borde=zeros(n+2,n+2);

%pilars coast

% r=2;
% borde=false(n+2,n+2);
% for i=2:4:n+1
%     borde=borde|((xx-i).^2+(yy-99).^2<r^2);
% end
% for i=2:5:n+1
%     borde=borde|((xx-i).^2+(yy-93).^2<r^2);
% end
% for i=2:6:n+1
%     borde=borde|((xx-i).^2+(yy-87).^2<r^2);
% end
% for i=2:8:n+1
%     borde=borde|((xx-i).^2+(yy-81).^2<r^2);
% end
% borde=1*borde;

%pilars separated from the coast

% r=2;
% borde=false(n+2,n+2);
% for i=2:5:n+1
%     borde=borde|((xx-i).^2+(yy-79).^2<r^2);
% end
% for i=2:6:n+1
%     borde=borde|((xx-i).^2+(yy-73).^2<r^2);
% end
% for i=2:8:n+1
%     borde=borde|((xx-i).^2+(yy-67).^2<r^2);
% end
% borde=1*borde;

Ah=zeros(n,n,5);
Au=zeros(n,n,5);
Av=zeros(n,n,5);

for i=1:n
    for j=1:n
        if borde(i+1,j+1)==0
            Ah(i,j,3)=1;
            Au(i,j,3)=1;
            Av(i,j,3)=1;
        else
            acum=0;
            if borde(i+1,j)==0
                Ah(i,j,1)=1;
                acum=acum+1;
            end
            if borde(i,j+1)==0
                Ah(i,j,2)=1;
                acum=acum+1;
            end
            if borde(i+2,j+1)==0
                Ah(i,j,4)=1;
                acum=acum+1;
            end
            if borde(i+1,j+2)==0
                Ah(i,j,5)=1;
                acum=acum+1;
            end
            if acum~=0
                Ah(i,j,1)=Ah(i,j,1)/acum;
                Ah(i,j,2)=Ah(i,j,2)/acum;
                Ah(i,j,4)=Ah(i,j,4)/acum;
                Ah(i,j,5)=Ah(i,j,5)/acum;
            else
                Ah(i,j,3)=1;
            end
            
            acum=0;
            if borde(i,j+1)==0
                Au(i,j,2)=-1;
                acum=acum+1;
            end
            if borde(i+2,j+1)==0
                Au(i,j,4)=-1;
                acum=acum+1;
            end
            if borde(i+1,j)==0
                Au(i,j,1)=1;
                acum=acum+1;
            end
            if borde(i+1,j+2)==0
                Au(i,j,5)=1;
                acum=acum+1;
            end
            if acum~=0
                Au(i,j,2)=Au(i,j,2)/acum;
                Au(i,j,4)=Au(i,j,4)/acum;
                Au(i,j,1)=Au(i,j,1)/acum;
                Au(i,j,5)=Au(i,j,5)/acum;
            else
                Au(i,j,3)=1;
            end
            
            acum=0;
            if borde(i+1,j)==0
                Av(i,j,1)=-1;
                acum=acum+1;
            end
            if borde(i+1,j+2)==0
                Av(i,j,5)=-1;
                acum=acum+1;
            end
            if borde(i,j+1)==0
                Av(i,j,2)=1;
                acum=acum+1;
            end
            if borde(i+2,j+1)==0
                Av(i,j,4)=1;
                acum=acum+1;
            end
            if acum~=0
                Av(i,j,1)=Av(i,j,1)/acum;
                Av(i,j,5)=Av(i,j,5)/acum;
                Av(i,j,2)=Av(i,j,2)/acum;
                Av(i,j,4)=Av(i,j,4)/acum;
            else
                Av(i,j,3)=1;
            end
        end
    end
end

% Initialize graphics

[surfplot,topography,boundary,top,start,stop] = initgraphics(n,dx);

% Results

t=[];
s=[];
M=[];

% Outer loop, restarts.

while get(stop,'value') == 0
    set(start,'value',0)
    
    H = ones(n+2,n+2)-Z;   U = zeros(n+2,n+2);  V = zeros(n+2,n+2);
    Hx = zeros(n+1,n+1); Ux = zeros(n+1,n+1); Vx = zeros(n+1,n+1);
    Hy = zeros(n+1,n+1); Uy = zeros(n+1,n+1); Vy = zeros(n+1,n+1);
    
    nstep = 0;
    
    % Inner loop, time steps.
    
    while get(start,'value')==0 && get(stop,'value')==0
        nstep = nstep + 1;
        
        % First half step
        
        % x direction
        i = 1:n+1;
        j = 1:n;
        
        % height
        Hx(i,j) = (H(i+1,j+1)+H(i,j+1))/2 - dt/(2*dx)*(U(i+1,j+1)-U(i,j+1));
        
        % x momentum
        Ux(i,j) = (U(i+1,j+1)+U(i,j+1))/2 -  ...
            dt/(2*dx)*((U(i+1,j+1).^2./H(i+1,j+1) + g/2*H(i+1,j+1).^2) - ...
            (U(i,j+1).^2./H(i,j+1) + g/2*H(i,j+1).^2))...
            -dt*nn^2*(U(i+1,j+1)+U(i,j+1)).*sqrt((U(i+1,j+1)+U(i,j+1)).^2+(V(i+1,j+1)+V(i,j+1)).^2)./(8*((H(i+1,j+1)+H(i,j+1))/2).^(10/3))...
            -dt*g*(H(i+1,j+1)+H(i,j+1)).*(Z(i+1,j+1)-Z(i,j+1))/(4*dx);
        
        % y momentum
        Vx(i,j) = (V(i+1,j+1)+V(i,j+1))/2 - ...
            dt/(2*dx)*((U(i+1,j+1).*V(i+1,j+1)./H(i+1,j+1)) - ...
            (U(i,j+1).*V(i,j+1)./H(i,j+1)))...
            -dt*nn^2*(V(i+1,j+1)+V(i,j+1)).*sqrt((U(i+1,j+1)+U(i,j+1)).^2+(V(i+1,j+1)+V(i,j+1)).^2)./(8*((H(i+1,j+1)+H(i,j+1))/2).^(10/3));
        
        % y direction
        i = 1:n;
        j = 1:n+1;
        
        % height
        Hy(i,j) = (H(i+1,j+1)+H(i+1,j))/2 - dt/(2*dy)*(V(i+1,j+1)-V(i+1,j));
        
        % x momentum
        Uy(i,j) = (U(i+1,j+1)+U(i+1,j))/2 - ...
            dt/(2*dy)*((V(i+1,j+1).*U(i+1,j+1)./H(i+1,j+1)) - ...
            (V(i+1,j).*U(i+1,j)./H(i+1,j)))...
            -dt*nn^2*(U(i+1,j+1)+U(i+1,j)).*sqrt((U(i+1,j+1)+U(i+1,j)).^2+(V(i+1,j+1)+V(i+1,j)).^2)./(8*((H(i+1,j+1)+H(i+1,j))/2).^(10/3));
        % y momentum
        Vy(i,j) = (V(i+1,j+1)+V(i+1,j))/2 - ...
            dt/(2*dy)*((V(i+1,j+1).^2./H(i+1,j+1) + g/2*H(i+1,j+1).^2) - ...
            (V(i+1,j).^2./H(i+1,j) + g/2*H(i+1,j).^2))...
            -dt*nn^2*(V(i+1,j+1)+V(i+1,j)).*sqrt((U(i+1,j+1)+U(i+1,j)).^2+(V(i+1,j+1)+V(i+1,j)).^2)./(8*((H(i+1,j+1)+H(i+1,j))/2).^(10/3))...
            -dt*g*(H(i+1,j+1)+H(i+1,j)).*(Z(i+1,j+1)-Z(i+1,j))/(4*dy);
        
        %threshhold
        Hx(Hx<0.1)=0.1;
        Hy(Hy<0.1)=0.1;
        
        % Second half step
        i = 2:n+1;
        j = 2:n+1;
        
        %friction terms
        Sfx=dt*nn^2*U(i,j).*sqrt(U(i,j).^2+V(i,j).^2)./H(i,j).^(10/3);
        Sfy=dt*nn^2*V(i,j).*sqrt(U(i,j).^2+V(i,j).^2)./H(i,j).^(10/3);
        
        %topography terms
        Sox=dt*g*H(i,j).*(Z(i+1,j)-Z(i-1,j))/(2*dx);
        Soy=dt*g*H(i,j).*(Z(i,j+1)-Z(i,j-1))/(2*dy);
        
        % height
        H(i,j) = H(i,j) - (dt/dx)*(Ux(i,j-1)-Ux(i-1,j-1)) - ...
            (dt/dy)*(Vy(i-1,j)-Vy(i-1,j-1));
        
        % x momentum
        U(i,j) = U(i,j) - (dt/dx)*((Ux(i,j-1).^2./Hx(i,j-1) + g/2*Hx(i,j-1).^2) - ...
            (Ux(i-1,j-1).^2./Hx(i-1,j-1) + g/2*Hx(i-1,j-1).^2)) ...
            - (dt/dy)*((Vy(i-1,j).*Uy(i-1,j)./Hy(i-1,j)) - ...
            (Vy(i-1,j-1).*Uy(i-1,j-1)./Hy(i-1,j-1)))...
            -Sfx-Sox;
        
        % y momentum
        V(i,j) = V(i,j) - (dt/dx)*((Ux(i,j-1).*Vx(i,j-1)./Hx(i,j-1)) - ...
            (Ux(i-1,j-1).*Vx(i-1,j-1)./Hx(i-1,j-1))) ...
            - (dt/dy)*((Vy(i-1,j).^2./Hy(i-1,j) + g/2*Hy(i-1,j).^2) - ...
            (Vy(i-1,j-1).^2./Hy(i-1,j-1) + g/2*Hy(i-1,j-1).^2))...
            -Sfy-Soy;
        
        %boundary conditions
        H(:,1) = 1+0.3*sin(2*pi*4*nstep*dt+2*pi*0.5*(1:n+2)*dx);      U(:,1) = U(:,2);       V(:,1) = V(:,2);
        H(:,n+2) = H(:,n+1);  U(:,n+2) = U(:,n+1);   V(:,n+2) = -V(:,n+1);
        H(1,:) = H(n+1,:);      U(1,:) = U(n+1,:);      V(1,:) = V(n+1,:);
        H(n+2,:) = H(2,:);  U(n+2,:) = U(2,:);  V(n+2,:) = V(2,:);
        
        H(2:n+1,2:n+1)=Ah(:,:,1).*(H(2:n+1,1:n)-1)+Ah(:,:,2).*(H(1:n,2:n+1)-1)+Ah(:,:,3).*(H(2:n+1,2:n+1)-1)+Ah(:,:,4).*(H(3:n+2,2:n+1)-1)+Ah(:,:,5).*(H(2:n+1,3:n+2)-1)+1;
        U(2:n+1,2:n+1)=Au(:,:,1).*U(2:n+1,1:n)+Au(:,:,2).*U(1:n,2:n+1)+Au(:,:,3).*U(2:n+1,2:n+1)+Au(:,:,4).*U(3:n+2,2:n+1)+Au(:,:,5).*U(2:n+1,3:n+2);
        V(2:n+1,2:n+1)=Av(:,:,1).*V(2:n+1,1:n)+Av(:,:,2).*V(1:n,2:n+1)+Av(:,:,3).*V(2:n+1,2:n+1)+Av(:,:,4).*V(3:n+2,2:n+1)+Av(:,:,5).*V(2:n+1,3:n+2);
        
        %threshhold
        H(H<0.1)=0.1;
        
        %Add result
        t=[t nstep*dt];
        s=[s std(H(:,92))*sqrt((n+2)*dx)];
        
        % Update plot
        if mod(nstep,nplotstep) == 0
            C = abs(U(i,j)) + abs(V(i,j));  % Color shows momemtum
            set(surfplot,'zdata',H(i,j)+Z(i,j),'cdata',C);
            luz=zeros(n,n,3);
            luz(:,:,1)=0.2*abs(U(i,j));
            luz(:,:,2)=0.2*abs(V(i,j));
            luz(:,:,3)=0.2*(abs(U(i,j))+abs(V(i,j)));
            set(topography,'zdata',Z(i,j),'cdata',luz);
            colorborde=zeros(n,n,3);
            colorborde(:,:,1)=0.3*borde(i,j);
            colorborde(:,:,2)=0.3*borde(i,j);
            colorborde(:,:,3)=0.3*borde(i,j);
            set(boundary,'zdata',2*borde(i,j),'cdata',colorborde);
            set(top,'string',sprintf('t = %6.2f',nstep*dt))
            drawnow
            M=[M getframe];
        end
        
        if all(all(isnan(H))), break, end  % Unstable, restart
    end
end
close(gcf)

end

% ------------------------------------

function [surfplot,topography,boundary,top,start,stop] = initgraphics(n,dx)
% INITGRAPHICS  Initialize graphics for waterwave.
% [surfplot,top,start,stop] = initgraphics(n)
% returns handles to a surface plot, its title, and two uicontrol toggles.

clf
shg
set(gcf,'numbertitle','off','name','Shallow_water','Position',[0 0 800 600])
x = (1:n)*dx;
surfplot = surf(x,x,ones(n,n),zeros(n,n));
axis off
hold on;
topography = surf(x,x,zeros(n,n),zeros(n,n));
boundary = surf(x,x,zeros(n,n),zeros(n,n));
grid off
axis([0 n*dx 0 n*dx -1 4])
caxis([-1 1])
shading interp
c = (1:64)'/64;
cyan = [0*c c.^2 c.*(2-c)];
colormap(cyan)
top = title('Click start');
start = uicontrol('position',[20 20 80 20],'style','toggle','string','start');
stop = uicontrol('position',[120 20 80 20],'style','toggle','string','stop');

end
