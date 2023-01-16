dt=0.01; %Time Step
T=20;
Time=[0:dt:T];


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


x=[1;1;1;1];
u=[0;0];


nn=size(A,1);
mm=size(B,2);
rr=size(C,1);


Ad=(A)*dt+eye(nn,nn);
Bd=B*dt;


Q = 2.3; 
R = 1; 
w = sqrt(Q)*randn(length(Time),1);
v = sqrt(R)*randn(length(Time),1);



for i=0:dt:T

x=Ad*x+Bd*u;
y=C*x;


X=[X x];
Y=[Y y];

end

plot(Time,X')

function cost=JJ(u)

for i=1:Np

    %x = xlim*dt*
    cost=cost+(x_estimate-y);
end
end