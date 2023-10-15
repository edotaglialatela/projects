clear
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% numerical solution for the response of a car's suspension system under a load
% the spring behaviour in the suspension system is elastic-plastic with hardening
% the default geometric values are fairly similar to a ferrari 296 gtb
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%1- ask for the parameters
prompt={'Rear spring constant [kg/cm = N/mm]',... %1
    'Front spring constant [kg/cm = N/mm]',... %2
    'Mass [kg]',... %3
    'Yield limit [MPa]',... %4
        'Rear damping factor [cP]',... %5
        'Front damping factor [cP]',... %6
        'Amplitude of force [N]',... %7
        'Pulsation of force factor [1/s]',... %8
        'Length 1 (front-center of mass) [m]',... %9
        'Length 2 (center of mass-applied force) [m]',... %10
        'Length 3 (applied force-rear end) [m]',... %11
        'Time interval [s]',... %12
        'Type of force: 1-constant, 2-alternating, 3-start',...%13
        'Initial position [m]',...%14
        'Initial velocity[m/s]',... %15
        'Time step'}; %16
name='Input dati';
numlines=1;
defaultanswer={'51','51','1470','1600','7000','7000','1960','1000','1.503','1.1','0.45','1','3','0','0','1e-4'};
options.Resize='on';
options.WindowStyle='normal';
options.Interpreter='tex';
answerf=inputdlg(prompt,name,numlines,defaultanswer,options);

M = str2double(cell2mat(answerf(3))); %[N/m]
k_1 = str2double(cell2mat(answerf(2)))*10e3; %[N/m] LADA 2101-2107 RIVA NOVA
k_1p= 0.33.*k_1; % k_1 hardening ratio
k_2 = str2double(cell2mat(answerf(1)))*10e3; %[N/m] (Koni TrackDay kit)
k_2p = 0.33.*k_2; % k_2 hardening ratio
fy = ((pi*(0.01^3))/(16*0.04*sqrt(3)*1.0625))*((str2double(cell2mat(answerf(4))))*(10^6));  %yielding force [N]
mu_1 = str2double(cell2mat(answerf(6))); %MED-6346 silicone gel 
mu_2 = str2double(cell2mat(answerf(5)));
%l_tot = 458cm, 260cm wheel to wheel
%58.5% weight on the rear, 41.5% weight on the front => c.m. approx. at
%271.7cm (59% of l_tot) (total weight = 1470kg)
l_1 = str2double(cell2mat(answerf(9)));
l_2 = str2double(cell2mat(answerf(10)));
l_3 = str2double(cell2mat(answerf(11)));
l_4 = l_1 + l_2;
F0 = str2double(cell2mat(answerf(7))); %[N]
J0 = (0.595*M*(1.05^2))+(0.405*M*(1.55.^2)); %[kg*m^2]
omega = str2double(cell2mat(answerf(8))); %[rad/s]
dof = 2 ;
endtime = str2double(cell2mat(answerf(12))); %[s]
period = 2*pi/omega; %[s] 
time_step=str2double(cell2mat(answerf(16)));

%2- setup matrices
W = [M,0;0,J0]; %mass matrix
eig_W = eig(W); %is W positive-definite?
K = [k_1+k_2,k_2*l_3-k_1*l_4;k_2*l_3-k_1*l_4,k_1*((l_4)^2)+k_2*((l_3)^2)]; %stiffness matrix
eig_K = eig(K); %is K positive-definite?
C = [mu_1+mu_2,mu_2*l_3-mu_1*l_4;mu_2*l_3-mu_1*l_4,mu_1*((l_4)^2)+mu_2*((l_3)^2)]; %damping matrix 
I = eye(dof,dof); %identity matrix
Z = zeros(dof,dof); %zero matrix
A = [W , C ; Z , I]; 
B = [ Z , K ; -I , Z];

%2- ODE well-posedness check
G_prime=-A\B;
spec_radius = max(abs(eig(G_prime)));
K0 = sqrt(spec_radius);

%4- numerically integrate EOMs
npassi=(endtime)/time_step;
initial_cond = [str2double(cell2mat(answerf(15))).*ones(dof,1) ; str2double(cell2mat(answerf(14))).*ones(dof,1)];
t=zeros(npassi,1);
Y=[str2double(cell2mat(answerf(15)));str2double(cell2mat(answerf(15)));str2double(cell2mat(answerf(14)));str2double(cell2mat(answerf(14)))];
YY=zeros(4,npassi);
dXP=zeros(2,npassi);
    for i=2:npassi %matlab matrices start with 1
        t(i)=t(i-1)+time_step;
        [FF(:,i)]=forcefunct(str2double(cell2mat(answerf(13))),F0,l_2,t(i),omega);
        FFF=FF(:,i);
        Y(:,i)=Y(:,i-1)+((A\(FFF-B*Y(:,i-1)))*time_step);
        YY(:,i)=Y(i);
        front_action(i)=(k_1*(YY(3,i)-((YY(4,i))*l_1)));
        back_action(i)=(k_2*(YY(3,i)+((YY(4,i))*(l_2+l_3))));
    if str2double(cell2mat(answerf(13)))==1 %constant
        if front_action(1,i)>=fy && back_action(1,i)>=fy
            if sign(front_action(1,i))==1 && sign(back_action(1,i))==1
                dXP(:,i)=[dXP(1,i-1)+((front_action(i)/k_1p)-(fy/k_1));dXP(2,i-1)+((back_action(i)/k_2p)-(fy/k_2))];
            elseif sign(front_action(1,i))==1 && sign(back_action(1))==-1
                dXP(:,i)=[dXP(1,i-1)+((front_action(i)/k_1p)-(fy/k_1));dXP(2,i-1)];
            elseif sign(front_action(1,i))==-1 && sign(back_action(1))==1
                dXP(:,i)=[dXP(1,i-1);dXP(2,i-1)+((back_action(i)/k_2p)-(fy/k_2))];
            elseif sign(front_action(1,i))==-1 && sign(back_action(1))==-1
                dXP(:,i)=dXP(:,i-1);
            end 
        elseif front_action(1,i)>=fy && back_action(1,i)<fy
            if sign(front_action(1,i))==1
                dXP(:,i)=[dXP(1,i-1)+((front_action(i)/k_1p)-(fy/k_1));dXP(1,i-1)];
            elseif sign(front_action(1,i))==-1 
                dXP(:,i)=dXP(:,i-1);
            end
        elseif front_action(1,i)<fy && back_action(1,i)>=fy
            if sign(back_action(1,i))==1
                dXP(:,i)=[dXP(1,i-1);dXP(2,i-1)+((back_action(i)/k_2p)-(fy/k_2))];
            elseif sign(back_action(1,i))==-1 
                dXP(:,i)=dXP(:,i-1);
            end
        end
    
    elseif str2double(cell2mat(answerf(13)))==2 | str2double(cell2mat(answerf(13)))==3 %alternating or start
        sgnF1(i)=sign(front_action(i)-front_action(i-1));
        sgnM2(i)=sign(back_action(i)-back_action(i-1));
        if front_action(i)>=fy && back_action(i)>=fy
            if sgnF1(i)==1 && sgnM2(i)==1
                dXP(:,i)=[dXP(1,i-1)+((front_action(i)/k_1p)-(fy/k_1));dXP(2,i-1)+((back_action(i)/k_2p)-(fy/k_2))];
            elseif sgnF1(i)==1 && sgnM2(i)==-1
                dXP(:,i)=[dXP(1,i-1)+((front_action(i)/k_1p)-(fy/k_1));dXP(2,i-1)];          
            elseif sgnF1(i)==-1 && sgnM2(i)==1
                dXP(:,i)=[dXP(1,i-1);dXP(2,i-1)+((back_action(i)/k_2p)-(fy/k_2))];
            elseif sgnF1(i)==-1 && sgnM2(i)==-1
                dXP(:,i)=dXP(:,i-1);
            end
        elseif front_action(i)>=fy && back_action(i)<fy
            if sgnF1(i)==1
                dXP(:,i)=[dXP(1,i-1)+((front_action(i)/k_1p)-(fy/k_1));dXP(2,i-1)];
            else
                dXP(:,i)=dXP(:,i-1);
            end
        elseif front_action(i)<fy && back_action(i)>=fy
            if sgnM2(i)==1
                dXP(:,i)=[dXP(1,i-1);dXP(2,i-1)+((back_action(i)/k_2p)-(fy/k_2))];
            else
                dXP(:,i)=dXP(:,i-1);
            end
        else
            dXP(:,i)=dXP(:,i-1);
        end
    end
    end

%5- plot results

plot(t,FF(1,:),"Color",[0.6350 0.0780 0.1840])
hold on
plot(t,FF(2,:),"Color",[0 0.4470 0.7410])
hold off
legend('F_1','F_2')
xlabel("t")
title("Actions")

plot(t,Y(1,:),"Color",[0.6350 0.0780 0.1840])
hold on
plot(t,Y(3,:),"Color",[0 0.4470 0.7410])
hold off
legend('v','x')
xlabel("t")
title("Response 1")
plot(t,Y(3,:),"Color",[0 0.4470 0.7410])
xlabel('t')
ylabel('x_{cm}(t)')

plot(t,Y(2,:),"Color",[0.6350 0.0780 0.1840])
hold on
plot(t,Y(4,:),"Color",[0 0.4470 0.7410])
hold off
legend('\omega','\theta')
xlabel("t")
title("Response 2")
plot(t,Y(4,:),"Color",[0 0.4470 0.7410])
xlabel('t')
ylabel('\vartheta(t)')

plot(Y(3,:),Y(1,:),"Color",[0.6350 0.0780 0.1840])
ylabel("v")
xlabel("x")
title("Phase space 1")

plot(Y(4,:),Y(2,:),"Color",[0.6350 0.0780 0.1840])
ylabel("\omega")
xlabel("\theta")
title("Phase space 2")

plot(t,dXP(1,:),"Color",[0.6350 0.0780 0.1840])
ylabel("x_1^p")
xlabel("t")
title("Plastic strain on the front springs")

plot(t,dXP(2,:),"Color",[0.6350 0.0780 0.1840])
ylabel("x_2^p")
xlabel("t")
title("Plastic strain on the back springs")