%%
% Please refer to Section 4, Simulation and discussion
% Xiangdong Yang, Liao Wu, Jinquan Li, Ken Chen, A minimal kinematic model for serial robot calibration using POE formula, Robotics and Computer-Integrated Manufacturing, Volume 30, Issue 3, June 2014, Pages 326-334
clc;
clear;

%% prepare data
s=0.001; %unit scale for distance, s=0.001 if use m as unit, s=1 if use mm as unit

%scara robot model, nominal parameters
w1=[0;0;1];     p1=[0;0;0]*s;
w2=[0;0;1];     p2=[250;0;0]*s;
w3=[0;0;0];     p3=[0;0;-1];
w4=[0;0;-1];    p4=[470;0;0]*s;
xi0=[twistCoord(w1,p1),twistCoord(w2,p2),twistCoord(w3,p3),twistCoord(w4,p4)];

%joint offset error
deq1=0;
deq2=0.02;
deq3=2*s;
deq4=0.02;

%actual joint twists
xi01=[0.0199900035972015;0;0.999800179914059;0;0.0130330000000000*1000*s;0;];
xi02=[0;0.000399999968000004;0.999999920000010;-0.000300000000000000*1000*s;-0.253990000161600*1000*s;0.000101596000064640*1000*s;];
xi03=[0;0;0;0.0199999568791395;0.0195999577415567;-0.999607844797830;];
xi04=[0.0407700195329210;0.0391700187663605;-0.998400478333784;-0.0266829837144079*1000*s;0.504558015646471*1000*s;0.0187056011887379*1000*s;];
xi00=[xi01,xi02,xi03,xi04];

N=10;%measurement number
vtheta=[rand(N,1)*2*pi,rand(N,1)*2*pi,rand(N,1)*1000*s,rand(N,1)*2*pi];%N groups of random joint positions

P00=[-100;-100;-100]*s;%point for validation
P01=[100;0;0]*s;%nominal position of point 1
P02=[0;100;0]*s;%nominal position of point 2
P03=[0;0;100]*s;%nominal position of point 3
PX=[P01,P02,P03];

gm=zeros(4,4,N);%measured end-effector poses
gn=zeros(4,4,N);%nominal end-effector poses

Pa1=zeros(3,N);%measured position of point 1
Pa2=zeros(3,N);%measured position of point 2
Pa3=zeros(3,N);%measured position of point 3

%simulate measurement data, add random error
for i=1:N
    Pa1(:,i)=[eye(3),zeros(3,1)]*fkSCARA(xi00,[vtheta(i,1)+deq1;vtheta(i,2)+deq2;vtheta(i,3)+deq3;vtheta(i,4)+deq4],4)*[P01;1]+(rand(3,1)*0.2-0.1)*s;
    Pa2(:,i)=[eye(3),zeros(3,1)]*fkSCARA(xi00,[vtheta(i,1)+deq1;vtheta(i,2)+deq2;vtheta(i,3)+deq3;vtheta(i,4)+deq4],4)*[P02;1]+(rand(3,1)*0.2-0.1)*s;
    Pa3(:,i)=[eye(3),zeros(3,1)]*fkSCARA(xi00,[vtheta(i,1)+deq1;vtheta(i,2)+deq2;vtheta(i,3)+deq3;vtheta(i,4)+deq4],4)*[P03;1]+(rand(3,1)*0.2-0.1)*s;
    PY=[Pa1(:,i),Pa2(:,i),Pa3(:,i)];
    
    [R,t,~,~,~]=Registration(PX,PY,eye(3),zeros(3,1),1);%use point based registration to get end-effector pose
    gm(1:3,1:3,i)=R;
    gm(1:3,4,i)=t;

    gn(:,:,i)=fkSCARA(xi00,[vtheta(i,1)+deq1;vtheta(i,2)+deq2;vtheta(i,3)+deq3;vtheta(i,4)+deq4],4);%nominal end-effector pose
end

%% calibration
[xiMinimal,dqMinimal,meanEMinimal,convergenceMinimal]=ScaraMinimal(xi0,vtheta,gm,10);%calibration with minimal model
[xiTraditional,dqTraditional,meanETraditional,convergenceTraditional]=ScaraTraditional(xi0,vtheta,gm,10);%calibration with traditional model

%% evaluation
Nt=50; %number of test data
error_before=zeros(Nt,1);%error before calibration
error_afterMinimal=zeros(Nt,1);%error after calibration with minimal model
error_afterTraditional=zeros(Nt,1);%error after calibration with traditional model

TestJointConfig=[rand(Nt,1)*2*pi,rand(Nt,1)*2*pi,rand(Nt,1)*1000*s,rand(Nt,1)*2*pi]; %Nt groups of test configurations

for i=1:Nt
    error_before(i)=norm((fkSCARA(xi0,TestJointConfig(i,:),4)-fkSCARA(xi00,TestJointConfig(i,:)+[deq1,deq2,deq3,deq4],4))*[P00;1])/s;
    error_afterMinimal(i)=norm((fkSCARA(xiMinimal,TestJointConfig(i,:)+dqMinimal',4)-fkSCARA(xi00,TestJointConfig(i,:)+[deq1,deq2,deq3,deq4],4))*[P00;1])/s;
    error_afterTraditional(i)=norm((fkSCARA(xiTraditional,TestJointConfig(i,:)+dqTraditional',4)-fkSCARA(xi00,TestJointConfig(i,:)+[deq1,deq2,deq3,deq4],4))*[P00;1])/s;
end

%% output result
mean(error_before)
mean(error_afterMinimal)
mean(error_afterTraditional)

max(error_before)
max(error_afterMinimal)
max(error_afterTraditional)