% A simple code for simulation of linear dynamic systems in State space format


%Time Step
dt=0.01; 
%Simulation final time
T=20;
% Time Vector
Time=[0:dt:T]; 


% Continous time state space representation of the system

A = [-0.0558   -0.9968    0.0802    0.0415
      0.5980   -0.1150   -0.0318         0
     -3.0500    0.3880   -0.4650         0
           0    0.0805    1.0000         0];

B = [ 0.0073         0
     -0.4750    0.0077
      0.1530    0.1430
           0         0];

C = [0     1     0     0
     0     0     0     1];

D = [0     0
     0     0];

% Initial condition of the system
x=[1;1;1;1];
u=[0;0];


% Number of states (nn), inputs (mm) and measurment sensors(rr)

nn=size(A,1);
mm=size(B,2);
rr=size(C,1);


% Discrete system with first order Euler discretization method
Ad=(A)*dt+eye(nn,nn);
Bd=B*dt;


%Place holder for saving states and measurments
X=[];
Y=[];

% Simulation loop
for i=0:dt:T

x=Ad*x+Bd*u;
y=C*x;


X=[X x];
Y=[Y y];

end

%Plot the states
plot(Time,X')
xlabel("Time")
ylabel("States")
% Plot the measurments
figure
plot(Time,Y')
xlabel("Measurments")
ylabel("Time")



function cost=JJ(u)

for i=1:Np

    %x = xlim*dt*
    cost=cost+(x_estimate-y);
end
end