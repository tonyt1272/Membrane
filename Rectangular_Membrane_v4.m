%%Anthony Thompson, 11/17/2018,  Comprehensive Exams Review
%Wave propagation on a rectangular membrane with fixed edges.
%Initial displacement is a guassian

%%Testing git version tracking%%

clear
clc
close all
close hidden

j=1i;
L_x=1;   
L_z=1;
sigma_x=L_x/20;%.05  %Width (standard dev) of the Gaussian initial disturbance
sigma_z=L_z/20;
mu_x=L_x/2;
mu_z=L_z/2;
h=.5;
Beta=.25; %Damping Coeffecient
T=70;    %Tension in the surface 
rho_S=1; % area density, (Kg/meter^2)
c=sqrt(T/rho_S); %wave speed on the membrane (meters/second)
prompt='Hit return when ready';
n_modes=19;  %number of modes in the simulation n=1 to n_modes, m=1 to n_modes
for m=1:n_modes
    for n=1:n_modes     %Generating mode frequencies
        w_mn(m,n)=c*pi*sqrt((n/L_x)^2+(m/L_z)^2);
    end
end

n=1;
m=1;
dx=.005;
dz=.005;
dt=.001;
x=0:dx:L_x;
z=0:dz:L_z;
[X,Z] = meshgrid(x,z);

%Initial = h*exp(-((X-mu_x).^2+(Z-mu_z).^2)/(2*sigma_x^2));  %The initial disturbance
Initial = h*exp(-(((X-mu_x).^2/(2*sigma_x^2))+((Z-mu_z).^2)/(2*sigma_z^2)));
%Initial = h*exp(-(((X-mu_x).^2)/(2*sigma_x^2)+(Z-mu_z).^2)/(2*sigma_x^2));  %The initial disturbance

for i=1:n_modes         %building the mode shapes 
    ModeX(:,:,i)=sin(i*pi*X/L_x);
    ModeZ(:,:,i)=sin(i*pi*Z/L_z);
end

Shape0=zeros(size(Initial));
figure(1);
hFig1 = figure(1);
set(hFig1, 'Position', [250 200 800 700])
xlim([-L_x,L_x])
ylim([-L_z,L_z])
t=0;
t_duration=.5;

for t=0:dt:t_duration
for n=1:n_modes
    for m=1:n_modes         %Applying amplitudes to modes and time progression
        A(n,m)=(4/(L_x*L_z))*sum(sum(ModeX(:,:,n).*ModeZ(:,:,m).*Initial)*dx*dz);
        Shape0(:,:)=Shape0(:,:)+A(n,m)*(ModeX(:,:,n).*ModeZ(:,:,m))*(cos(w_mn(n,m)*t))*exp(-Beta*t);

    end
end
    
    surf(X,Z,Shape0)
    zlim([-h,h])
    pause(.01)
    if t==0
        title(['Thin Membrane,Fixed Edge,Damped, modes=',num2str(n_modes),])
        xlabel('X')
        ylabel('Z')
        str = input(prompt,'s');
    end
    Shape0=zeros(size(Initial));
end


figure(2)
hFig2 = figure(2);
set(hFig2, 'Position', [450 200 800 700])
surf(X,Z,Initial)
title('The Initial Gaussian Displacement')
zlim([-1.1*h,1.1*h])
xlabel('X')
ylabel('Z')



