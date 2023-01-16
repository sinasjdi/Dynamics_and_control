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



%{
A=[0	1	0	0;...
0	-0.100000000000000	3	0;...
0	0	0	1;...
0	-0.500000000000000	30	0];

B=[0;2;0;5];

C=[1	0	0	0;...
0	0	1	0];


x=[1;1;1;1];
u=[0];
%}
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

mu_model_uncertainty=0.00;
sigma_model_uncertainty=0.001;


mu_measurment_uncertainty=0;
sigma_measurment_uncertainty=0.1;

Q = 0.5
R = 1;

[Kd,S,e] = dlqr(Ad,Bd,Q,R);

% Simulation loop
for i=0:dt:T

system_uncrtainty=normrnd(mu_model_uncertainty,sigma_model_uncertainty,size(mm,1));
measurement_noise=normrnd(mu_measurment_uncertainty,sigma_measurment_uncertainty,size(rr,1));
u=-Kd*x
x=Ad*x+Bd*u+system_uncrtainty;
y=C*x+measurement_noise;



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



