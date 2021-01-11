function [ xi,dq,meanE,convergence ] = ScaraMinimal( xi0,vtheta,gm,M )
%[xi,dq,meanE,convergence] = ScaraMinimal(xi0,vtheta,gm,M) : Minimal calibration model for scara robot
%
%Input:
%   xi0: nominal twists; 
%   vtheta: joint positions; 
%   gm: Actual measured end-effector poses
%   M: iteration steps
%
%Output:
%   xi: calibrated twists; 
%   dq: calibrated joint offsets; 
%   meanE: mean error
%   convergence: [Residual error, Variable increment]
%
%Please refer to 'Xiangdong Yang, Liao Wu, Jinquan Li, Ken Chen, A minimal
%kinematic model for serial robot calibration using POE formula, Robotics and Computer-Integrated Manufacturing, Volume 30, Issue 3, June 2014, Pages 326-334'

%memory allocation
    xi=xi0;
    dq=zeros(4,1);
    N=size(vtheta,1);

    gn=zeros(4,4,N);%nominal end-effector poses
    dg=zeros(4,4,N);%pose error
    vLog=zeros(6,N);%pose error in twist form
    for i=1:N
        gn(:,:,i)=fkSCARA(xi,vtheta(i,:),4);
        dg(:,:,i)=gm(:,:,i)/gn(:,:,i);
        vLog(:,i)=vlog(dg(:,:,i));
    end
    error=zeros(3,1);
    for i=1:N
        error=error+[norm(vLog(:,i));norm(vLog(4:6,i));norm(vLog(1:3,i))];%[total error;position error;rotation error]
    end
    meanE(:,1)=error/N;

    convergence=zeros(M,2);
    for m=1:M
        for i=1:N
            gn(:,:,i)=fkSCARA(xi,vtheta(i,:),4);
            dg(:,:,i)=gm(:,:,i)/gn(:,:,i);
            vLog(:,i)=vlog(dg(:,:,i));
        end
        simY=zeros(6*N,1);
        for i=1:N
            simY(6*i-5:6*i,1)=vLog(:,i);
        end
        %transformation from base frame to link frame
        t1=twistframe(xi(:,1));
        t2=twistframe(xi(:,2));
        t3=twistframe(xi(:,3));
        t4=twistframe(xi(:,4));
        
        %Jacobian matrix
        for k=0:N-1
            simJ(1+6*k:6+6*k,1:7)=dexp(xi(:,1),vtheta(k+1,1));
            simJ(1+6*k:6+6*k,8:14)=adM(se3Exp(vtheta(k+1,1)*xi(:,1)))*dexp(xi(:,2),vtheta(k+1,2));
            simJ(1+6*k:6+6*k,15:21)=adM(se3Exp(vtheta(k+1,1)*xi(:,1))*se3Exp(vtheta(k+1,2)*xi(:,2)))*dexp(xi(:,3),vtheta(k+1,3));
            simJ(1+6*k:6+6*k,22:28)=adM(se3Exp(vtheta(k+1,1)*xi(:,1))*se3Exp(vtheta(k+1,2)*xi(:,2))*se3Exp(vtheta(k+1,3)*xi(:,3)))*dexp(xi(:,4),vtheta(k+1,4));
        end

        Brot=[1, 0, 0, 0;
            0, 1, 0, 0;
            0, 0, 0, 0;
            0, 0, 0, 1;
            0, 0, -1, 0;
            0, 0, 0, 0];
        Btrl=[0,0;
            0, 0;
            0, 0;
            1, 0;
            0, 1;
            0, 0];
        simJ1=zeros(6*N,18);
        simJ1(:,1:4)=simJ(:,1:6)*adM(t1)*Brot;
        simJ1(:,5)=simJ(:,7);
        simJ1(:,6:9)=simJ(:,8:13)*adM(t2)*Brot;
        simJ1(:,10)=simJ(:,14);
        simJ1(:,11:12)=simJ(:,15:20)*adM(t3)*Btrl;
        simJ1(:,13)=simJ(:,21);
        simJ1(:,14:17)=simJ(:,22:27)*adM(t4)*Brot;
        simJ1(:,18)=simJ(:,28);

        dp=simJ1\simY;
        
        %update
        omega1=[dp(1);dp(2);(1-dp(1)*dp(1)-dp(2)*dp(2))^0.5];
        q1=[dp(3);dp(4);0];
        v1=cross(q1,omega1);
        xi(:,1)=adM(t1)*[omega1;v1];
        omega2=[dp(6);dp(7);(1-dp(6)*dp(6)-dp(7)*dp(7))^0.5];
        q2=[dp(8);dp(9);0];
        v2=cross(q2,omega2);
        xi(:,2)=adM(t2)*[omega2;v2];
        v3=[dp(11);dp(12);(1-dp(11)*dp(11)-dp(12)*dp(12))^0.5];
        xi(:,3)=adM(t3)*[0;0;0;v3];
        omega4=[dp(14);dp(15);(1-dp(14)*dp(14)-dp(15)*dp(15))^0.5];
        q4=[dp(16);dp(17);0];
        v4=cross(q4,omega4);
        xi(:,4)=adM(t4)*[omega4;v4];

        vtheta(:,1)=vtheta(:,1)+dp(5);
        vtheta(:,2)=vtheta(:,2)+dp(10);
        vtheta(:,3)=vtheta(:,3)+dp(13);
        vtheta(:,4)=vtheta(:,4)+dp(18);

        dq=dq+dp([5,10,13,18]);
        
        for i=1:N
            gn(:,:,i)=fkSCARA(xi,vtheta(i,:),4);
            dg(:,:,i)=gm(:,:,i)/gn(:,:,i);
            vLog(:,i)=vlog(dg(:,:,i));
        end
        error=zeros(3,1);
        for i=1:N
            error=error+[norm(vLog(:,i));norm(vLog(4:6,i));norm(vLog(1:3,i))];
        end
        meanE(:,m+1)=error/N;
        convergence(m,1)=norm(simY); %residual error
        convergence(m,2)=norm(dp); %variable increment
    end
end
