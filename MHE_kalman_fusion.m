number_of_sensors=5;
number_of_states=6
global N_e R_matrix meas_horizon C
RR=0.64
QQ=100000

RR_1=0.5
RR_2=0.5
RR_3=0.5
RR_4=0.5
RR_5=0.5

QQ = QQ*eye(number_of_states,number_of_states);


R_1=RR_1*eye(number_of_states,number_of_states);


R_2=RR_2*eye(number_of_states,number_of_states);
R_3=RR_3*eye(number_of_states,number_of_states);
R_4=RR_4*eye(number_of_states,number_of_states);
R_5=RR_5*eye(number_of_states,number_of_states);

%R_sembel=[R_1,R_2,R_3,R_4,R_5]'
% [tag1 state 1; tag1 state2; tag1 state3;....tag 2 state1....tag6 state 6]
R_matrix=diag([RR_1,RR_1,RR_1,RR_1,RR_1,RR_1,RR_2,RR_2,RR_2,RR_2,RR_2,RR_2,RR_3,RR_3,RR_3,RR_3,RR_3,RR_3,RR_4,RR_4,RR_4,RR_4,RR_4,RR_4,RR_5,RR_5,RR_5,RR_5,RR_5,RR_5])
%R_matrix=diag([0.5,RR_1,RR_1,RR_1,RR_1,RR_1,10,RR_2,RR_2,RR_2,RR_2,RR_2,4,RR_3,RR_3,RR_3,RR_3,RR_3,RR_4,RR_4,RR_4,RR_4,RR_4,RR_4,RR_5,RR_5,RR_5,RR_5,RR_5,RR_5])
C1=eye(number_of_states,number_of_states);
C=repmat(C1,number_of_sensors,1);

A=eye(number_of_states);

N=1000;

x1 =  sin(5*linspace(0, 10, N)*pi)';
x2 =  sin(5*linspace(0, 10, N)*pi)';
x3 =  sin(5*linspace(0, 10, N)*pi)';
x4 =  sin(5*linspace(0, 10, N)*pi)';
x5 =  sin(5*linspace(0, 10, N)*pi)';
x6 =  sin(5*linspace(0, 10, N)*pi)';

x=[x1, x2,x3,x4,x5,x6]';
size(x)



noise_1_1=normrnd(1,RR_1,[1,N]);
noise_2_1=normrnd(2,RR_2,[1,N]);
noise_3_1=normrnd(20,RR_3,[1,N]);
noise_4_1=normrnd(20,RR_4,[1,N]);
noise_5_1=normrnd(20,RR_5,[1,N]);
noise_6_1=normrnd(20,RR_5,[1,N]);

noise_tag_1=[noise_1_1;noise_2_1;noise_3_1;noise_4_1;noise_5_1;noise_6_1];


noise_1_2=normrnd(1,RR_1,[1,N]);
noise_2_2=normrnd(2,RR_2,[1,N]);
noise_3_2=normrnd(20,RR_3,[1,N]);
noise_4_2=normrnd(20,RR_4,[1,N]);
noise_5_2=normrnd(20,RR_5,[1,N]);
noise_6_2=normrnd(20,RR_5,[1,N]);

noise_tag_2=[noise_1_2;noise_2_2;noise_3_2;noise_4_2;noise_5_2;noise_6_2];

noise_1_3=normrnd(1,RR_1,[1,N]);
noise_2_3=normrnd(2,RR_2,[1,N]);
noise_3_3=normrnd(20,RR_3,[1,N]);
noise_4_3=normrnd(20,RR_4,[1,N]);
noise_5_3=normrnd(20,RR_5,[1,N]);
noise_6_3=normrnd(20,RR_5,[1,N]);

noise_tag_3=[noise_1_3;noise_2_3;noise_3_3;noise_4_3;noise_5_3;noise_6_3];


noise_1_4=normrnd(1,RR_1,[1,N]);
noise_2_4=normrnd(2,RR_2,[1,N]);
noise_3_4=normrnd(20,RR_3,[1,N]);
noise_4_4=normrnd(20,RR_4,[1,N]);
noise_5_4=normrnd(20,RR_5,[1,N]);
noise_6_4=normrnd(20,RR_5,[1,N]);

noise_tag_4=[noise_1_4;noise_2_4;noise_3_4;noise_4_4;noise_5_4;noise_6_4];



noise_1_5=normrnd(1,RR_1,[1,N]);
noise_2_5=normrnd(2,RR_2,[1,N]);
noise_3_5=normrnd(20,RR_3,[1,N]);
noise_4_5=normrnd(20,RR_4,[1,N]);
noise_5_5=normrnd(20,RR_5,[1,N]);
noise_6_5=normrnd(20,RR_5,[1,N]);

noise_tag_5=[noise_1_5;noise_2_5;noise_3_5;noise_4_5;noise_5_5;noise_6_5];












noise=[noise_tag_1;noise_tag_2;noise_tag_3;noise_tag_4;noise_tag_5];



xx=repmat(x,[number_of_sensors,1]);

measurements=xx+noise;

plot(measurements');


P=1*eye(6);
x_est=x(:,1);
x_mhe=x(:,1);
PP=[];

X_EST=[];
X_MHE=[];
N_e=4;
for i=1:N
z=measurements(:,i);
[x_est,P]=kalman_fuser(x_est,A,C,P,QQ,R_matrix,z);

if(i>N_e)
meas_horizon=measurements(:,i-N_e:i);
x_mhe=mhe(x_mhe);
X_MHE=[X_MHE x_mhe];
end    


X_EST=[X_EST x_est];
PP(:,:,i)=P;

end

plot_kalman_results(X_EST,measurements)


function x_mhe=mhe(x_mhe)



[x_mhe,fval] = fmincon(@JJ,x_mhe,[],[]);

end

function [x_est,P]=kalman_fuser(x_est,A,C,P,QQ,R_matrix,z)
x_est = A* x_est;
P = A * P * A' + QQ;
y=z-C*x_est;
K = P * C' / (C * P * C' + R_matrix);
x_est = x_est + K * y;
P = (eye(size(A)) - K * C) * P;

end


function plot_kalman_results(X_EST,measurements)
plot(X_EST(:,:)')
figure
hold on
plot(measurements(1,:))
plot(measurements(2,:))
plot(measurements(3,:))
plot(measurements(4,:))
plot(measurements(5,:))
end


function cost=JJ(x_mhe)
global N_e R_matrix meas_horizon C
cost=0;
for i=1:N_e

    %x = xlim*dt*
    measure_temp=meas_horizon(:,i);
    cost=cost+(C*x_mhe-measure_temp)'*1*(C*x_mhe-measure_temp);
end
end