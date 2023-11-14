function[time,timedata]=DETM
clc;
close all;
clear;
h=0.001;
t=0:h:10; %时间
A1=[-0.1 50;-1 -10];B1=[0;1];C1=[0.3 0.5];D1=0;L1=[1 0];
A2=[-4.6 50;-1 -10];B2=[0;1];C2=[0.6 0.3];D2=0;L2=[1 0];
Nor=@(S,t)t'*S*t;
Af1=[-11.65,1.38;7.55,-7.99];Af2=[-27.13,2.31;-1.21,-6.15]; Bf1=[-0.06;0.01];Bf2=[-0.01;0.07];
Cf1=[5.85,-4.31];Cf2=[27.22,-4.11];
Df1=0;Df2=0;
O=[0.32];
T=0.2;
theta=4; lambda=1; r1=2;r2=0.9;
for k=1:length(t)
    if k==1
        tau(k)=1;
        x(:,k)=[2;3];
        y(:,k)=L1*x(:,k); 
%        y(:,k)=[2;2];
        xf(:,k)=[0.1;0.2];
        eta(k)=0.1;
        v(k)=0.1;
        w(k)=0;
        Num=1;
        time(Num)=k;
        j=1;
    else
            cc= x(1,k-1);
        if (-3 < cc & cc< 0)
            mem1=(cc+3)/3;
            elseif (0 <cc & cc< 2)
            mem1=(3-cc)/3;
            else
            mem1=0;
            end
            mem2=1-mem1;
            w(k)=5*exp(-1*t(k))*sin(5*t(k));                    
            x(1,k)=x(1,k-1)+h*(mem1*(A1(1,1)*x(1,k-1)+A1(1,2)*x(2,k-1))+mem2*(A2(1,1)*x(1,k-1)+A2(1,2)*x(2,k-1))+mem1*B1(1,1)*w(k-1)+mem2*B2(1,1)*w(k-1));
            x(2,k)=x(2,k-1)+h*(mem1*(A1(2,1)*x(1,k-1)+A1(2,2)*x(2,k-1))+mem2*(A2(2,1)*x(1,k-1)+A2(2,2)*x(2,k-1))+mem1*B1(2,1)*w(k-1)+mem2*B2(2,1)*w(k-1));                     
            y(:,k)=mem1*(L1(1,1)*x(1,k)+L1(1,2)*x(2,k))+mem2*(L2(1,1)*x(1,k)+L2(1,2)*x(2,k));      
            z(:,k)=mem1*(C1(1,1)*x(1,k)+C1(1,2)*x(2,k))+mem2*(C2(1,1)*x(1,k)+C2(1,2)*x(2,k));                           
            r(k)=y(:,k)-y(:,k-tau(k-1));
            R(k)=(r(k)'* r(k))/0.32;           
              if R(k)<1
                eta(k)=(exp(-r1*R(k))+1)*eta(k-1);
            elseif R(k)==1
                eta(k)=eta(k-1);
            elseif R(k)>1
                    eta(k)=exp(-r2*(R(k)-1))*eta(k-1);
              end            
            if eta(k)>0.5
                 eta(k)=0.5;
             elseif eta(k)<0.05
                 eta(k)=0.05;
            end       
        if tau(k-1)<j*T/h 
            v(k)=v(k-1)+h*(-lambda*v(k-1)+eta(k-1)*Nor(O,y(:,k-tau(k-1)+(j-1)*T/h ))-Nor(O,y(:,k-tau(k-1)+(j-1)*T/h)-y(:,k-tau(k-1))));
          eta(k)=eta(k-1);
         tau(k)=tau(k-1)+1;         
        else 
            if v(k-1)+theta*(eta(k)*(y(:,k)'*O*y(:,k))-(y(:,k)-y(:,k-tau(k-1)))'*O*(y(:,k)-y(:,k-tau(k-1))))>=0
            v(k)=v(k-1)+h*(-lambda*v(k-1)+eta(k-1)*Nor(O,y(:,k))-Nor(O,y(:,k )-y(:,k-tau(k-1))));
            tau(k)=tau(k-1)+1;  j=j+1;         
        else  
            j=1;
          tau(k)=1;   Num=Num+1; 
          time(Num)=k; 
          timedata(Num)=tau(k-1); 
          v(k)=v(k-1)+h*(-1*v(k-1)+eta(k-1)*Nor(O,y(:,k-1)));
            end
        end    
        xf(1,k)=xf(1,k-1)+h*(mem1*(Af1(1,1)*xf(1,k-1)+Af1(1,2)*xf(2,k-1))+mem2*(Af2(1,1)*xf(1,k-1)+Af2(1,2)*xf(2,k-1))+mem1*Bf1(1,1)*y(1,k-tau(k-1))+mem2*Bf2(1,1)*y(1,k-tau(k-1)));
        xf(2,k)=xf(2,k-1)+h*(mem1*(Af1(2,1)*xf(1,k-1)+Af1(2,2)*xf(2,k-1))+mem2*(Af2(2,1)*xf(1,k-1)+Af2(2,2)*xf(2,k-1))+mem1*Bf1(2,1)*y(1,k-tau(k-1))+mem2*Bf2(2,1)*y(1,k-tau(k-1)));
        zf(:,k)=mem1*(Cf1(1,1)*xf(1,k)+Cf1(1,2)*xf(2,k))+mem2*(Cf2(1,1)*xf(1,k)+Cf2(1,2)*xf(2,k));    
    end
end

figure (1);
stem(time*h,timedata*h);
axis([0 10 0 1.8]);
xlabel('t'); 
ylabel('The release time interval');

figure (2);
plot(t,z,'k',t,zf,'r','LineWidth',2);
xlabel('t'); 
ylabel('z(t) and z_f(t)');
legend({'$z(t)$','$z_f(t)$'},'interpreter','latex', 'Fontsize',12);
set(1,'position',[100 300 344*2 179*1.65]);

figure (3);
subplot(1,2,1);
plot(t,eta,'LineWidth',2)
legend({'$\eta$'},'interpreter','latex', 'Fontsize',14); 
xlabel('t'); 
ylabel('The trajectory of trigger threshold $\eta\left(t_{k} h\right)$','Interpreter','latex');
 axis([0 10 0 0.55]);
 grid on;
subplot(1,2,2);
plot(t,v,'LineWidth',2)
legend({'$\nu$'},'interpreter','latex', 'Fontsize',14); 
xlabel('t'); 
ylabel('The trajectory of dynamic variable $\nu$','Interpreter','latex');




