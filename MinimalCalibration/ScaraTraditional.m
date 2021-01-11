function [ xi,dq,meanE,convergence ] = ScaraTraditional( xi0,vtheta,gm,M)
%[xi,dq,meanE,convergence] = ScaraTraditional(xi0,vtheta,gm,M) : Traditional calibration model for scara robot
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

    xi=xi0;
    dq=zeros(4,1);
    N=size(vtheta,1);

    gn=zeros(4,4,N);
    dg=zeros(4,4,N);
    vLog=zeros(6,N);
    for i=1:N
        gn(:,:,i)=fkSCARA(xi,vtheta(i,:),4);
        dg(:,:,i)=gm(:,:,i)/gn(:,:,i);
        vLog(:,i)=vlog(dg(:,:,i));
    end
    error=zeros(3,1);
    for i=1:N
        error=error+[norm(vLog(:,i));norm(vLog(4:6,i));norm(vLog(1:3,i))];
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
        for k=0:N-1
            simJ(1+6*k:6+6*k,1:7)=dexp(xi(:,1),vtheta(k+1,1));
            simJ(1+6*k:6+6*k,8:14)=adM(se3Exp(vtheta(k+1,1)*xi(:,1)))*dexp(xi(:,2),vtheta(k+1,2));
            tempSIMJ=adM(se3Exp(vtheta(k+1,1)*xi(:,1))*se3Exp(vtheta(k+1,2)*xi(:,2)))*dexp(xi(:,3),vtheta(k+1,3));
            simJ(1+6*k:6+6*k,15:18)=tempSIMJ(:,4:7);
            simJ(1+6*k:6+6*k,19:25)=adM(se3Exp(vtheta(k+1,1)*xi(:,1))*se3Exp(vtheta(k+1,2)*xi(:,2))*se3Exp(vtheta(k+1,3)*xi(:,3)))*dexp(xi(:,4),vtheta(k+1,4));
        end
        dp=simJ\simY;

        xi(:,1)=xi(:,1)+dp(1:6);
        xi(1:3,1)=xi(1:3,1)/norm(xi(1:3,1));
        xi(4:6,1)=xi(4:6,1)-xi(1:3,1)'*xi(4:6,1)/(xi(1:3,1)'*xi(1:3,1))*xi(1:3,1);

        xi(:,2)=xi(:,2)+dp(8:13);
        xi(1:3,2)=xi(1:3,2)/norm(xi(1:3,2));
        xi(4:6,2)=xi(4:6,2)-xi(1:3,2)'*xi(4:6,2)/(xi(1:3,2)'*xi(1:3,2))*xi(1:3,2);

        xi(4:6,3)=xi(4:6,3)+dp(15:17);
        xi(4:6,3)=xi(4:6,3)/norm(xi(4:6,3));

        xi(:,4)=xi(:,4)+dp(19:24);
        xi(1:3,4)=xi(1:3,4)/norm(xi(1:3,4));
        xi(4:6,4)=xi(4:6,4)-xi(1:3,4)'*xi(4:6,4)/(xi(1:3,4)'*xi(1:3,4))*xi(1:3,4);

        vtheta(:,1)=vtheta(:,1)+dp(7);
        vtheta(:,2)=vtheta(:,2)+dp(14);
        vtheta(:,3)=vtheta(:,3)+dp(18);
        vtheta(:,4)=vtheta(:,4)+dp(25);

        dq=dq+dp([7,14,18,25]);
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

