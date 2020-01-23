%---------------------------------------------------------------
%    use the explicit scheme to solve the equation
%    u_t(x,t) = (a(u)u)_x(x,t),          xl = 0 < x < xr = 1, 0 < t < tf = 0.1
%    u(x,0) = f(x) = sin(2xpi),                   
%    u(0,t) = u(1,t) = 0.
%
%    Where we assume that there exists constants a_1 and a_2 s.t. 
%    0 < a_1 < a(u) < a_2
%    
%    and let a(u) be defined by,
%    a(u) = (1+2u^2)/(1+u^2).
%  
%    then we solve this program by the explicit scheme:
%    u(n,m+1) = u(n,m)+mu*(((2*u(n,m)*(u(n+1,m)-u(n,m))^2)/(1+(u(n,m)^2))^2)+((1+2*u(n,m)^2)*(u(n+1,m)-2*u(n,m)+u(n-1,m))/(1+u(n,m)^2)))
%---------------------------------------------------------------
clear all;                 % clear all variables in memory
close all;

xl=0; xr=1;                 % x domain [xl,xr]
N = 50;                     % N: number of mesh points for x
dx = (xr-xl) /(N-1);        % dx: mesh size
dt=0.0001;    				% dt: time step size
tmax = 0.1;                 % final simulation time
mu = dt/(dx)^2;				
Nt=round(tmax/dt);          % Nt: number of time steps

% Evaluate the initial conditions
x = xl : dx : xr;       	% generate the grid point       
f = sin(pi*2*x);  

% store the solution at all grid points for all time steps
u = zeros(N,Nt);   
u(:,1)=f(:);

% Find the approximate solution at each time step
for m = 1:Nt
    t = m*dt; 
    % boundary condition at left side        
    gl = 0;   
    % boundary condition at right side
    gr = 0;  
    for n=2:N-1    
      u(n,m+1)=u(n,m)+mu*(((2*u(n,m)*(u(n+1,m)-u(n,m))^2)/(1+(u(n,m)^2))^2)+((1+2*u(n,m)^2)*(u(n+1,m)-2*u(n,m)+u(n-1,m))/(1+u(n,m)^2)));
    end
       u(1,m+1) = gl;   
       u(N,m+1) = gr; 
end

% Plotting the result
tt = 0 : dt : Nt*dt;
figure(1)
surf(x,tt, u','EdgeColor','none');     
camlight left;
lighting phong
xlabel('x')
ylabel('t')
zlabel('u')
title('Numerical solution of 1-D nonlinear heat equation')
rotate3d on

% Animation
m1=VideoWriter('nonlinear.avi');
open(m1);
for m=1: round(Nt/300): Nt
    uu=u(:,m)';
    figure(2);
    plot(x,uu);
    xlim([0,0.5]);
    ylim([0,1]);
    xlabel('x')
    ylabel('u')
    t=m*dt;
    title(sprintf('t=%5.2f',t))
    frame=getframe;
    writeVideo(m1,frame);
    drawnow
    pause(0.1);
end
close(m1)
