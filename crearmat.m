function []=crearmat()
%

n=64;                  % grid size
g=9.8;                 % gravitational constant
dt=0.01;               % hardwired timestep
dx=1;
dy=1;
nplotstep=10;           % plot interval
nn=0.5;

x=1:n+2;
y=1:n+2;
[yy,xx]=meshgrid(y,x);
Z=yy/(10*(n+2));

r=2;

%borde=zeros(n+2,n+2);

borde=1*((((xx-33).^2+(yy-33).^2)<r^2)|(((xx-17).^2+(yy-39).^2)<r^2)|(((xx-49).^2+(yy-39).^2)<r^2)|(((xx-12).^2+(yy-45).^2)<r^2)|(((xx-33).^2+(yy-45).^2)<r^2)|(((xx-54).^2+(yy-45).^2)<r^2)|(((xx-8).^2+(yy-51).^2)<r^2)|(((xx-24).^2+(yy-51).^2)<r^2)|(((xx-40).^2+(yy-51).^2)<r^2)|(((xx-56).^2+(yy-51).^2)<r^2)|(((xx-7).^2+(yy-57).^2)<r^2)|(((xx-20).^2+(yy-57).^2)<r^2)|(((xx-33).^2+(yy-57).^2)<r^2)|(((xx-46).^2+(yy-57).^2)<r^2)|(((xx-59).^2+(yy-57).^2)<r^2));

clear x y xx yy r;
save parametros;

end