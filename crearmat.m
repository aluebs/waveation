function []=crearmat()
%

n=100;                  % grid size
g=9.8;                 % gravitational constant
dt=0.001;               % hardwired timestep
dx=0.1;
dy=0.1;
nplotstep=10;           % plot interval
nn=0.5;

x=1:n+2;
y=1:n+2;
[yy,xx]=meshgrid(y,x);

%Z=yy/(10*(n+2));

Z=0.9*(yy-(n-8))/10;
Z(Z<0)=0;

r=2;

%borde=zeros(n+2,n+2);

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

borde=false(n+2,n+2);
for i=2:5:n+1
    borde=borde|((xx-i).^2+(yy-79).^2<r^2);
end
for i=2:6:n+1
    borde=borde|((xx-i).^2+(yy-73).^2<r^2);
end
for i=2:8:n+1
    borde=borde|((xx-i).^2+(yy-67).^2<r^2);
end
borde=1*borde;

clear x y xx yy r;
save parametros;

end