clear;
clc;
close all;
%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%
dt=0.001;   %
t=10/dt+1;  %
A1=[-0.1 50;-1 -10]; B1=[0;1];C1=[0.3 0.5];D1=0;L1=[1 0]; 
A2=[-4.6 50;-1 -10];B2=[0;1];C2=[0.6 0.3];D2=0;L2=[1 0];I=1;
g=0.5; %最大阈值
M=0.2; %% 拉姆达
B=0.5; 
C=1; 
tauM=0.1; 
setlmis([]);
%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%
P1 = lmivar(1,[2,1]);
Q = lmivar(1,[2,1]);
W  = lmivar(1,[2,1]);

Af11= lmivar(1,[2,1]);
Bf11= lmivar(2,[2,1]);
Cf11= lmivar(2,[1,2]);

Af22= lmivar(1,[2,1]);
Bf22= lmivar(2,[2,1]);
Cf22= lmivar(2,[1,2]);

% Df1= lmivar(1,[1,1]);
R  = lmivar(1,[2,1]);
U  = lmivar(1,[2,1]);
O  = lmivar(1,[1,1]);
x = lmivar(1,[1,1]);

lmiterm([1 1 1 P1],1,A1,'s');
lmiterm([1 1 1 Q],1,1);
lmiterm([1 1 1 R],-1,1);

lmiterm([1 2 1 W],1,A1);
lmiterm([1 2 1 -Af11],1,1);

lmiterm([1 2 2 -Af11],1,1);
lmiterm([1 2 2 Af11],1,1);

lmiterm([1 3 1 -Bf11],L1',1);
lmiterm([1 3 1 R],1,1);
lmiterm([1 3 1 U],-1,1);

lmiterm([1 3 2 -Bf11],L1',1);

lmiterm([1 3 3 R],-2,1);
lmiterm([1 3 3 U],1,1);
lmiterm([1 3 3 -U],1,1);
lmiterm([1 3 3 O],g*L1',L1);

lmiterm([1 4 1 U],1,1);

lmiterm([1 4 3 R],1,1);
lmiterm([1 4 3 U],-1,1);

lmiterm([1 4 4 Q],-1,1);
lmiterm([1 4 4 R],-1,1);

lmiterm([1 5 1 -Bf11],-1,1);

lmiterm([1 5 2 -Bf11],1,1);

lmiterm([1 5 5 O],-1,1);

lmiterm([1 6 1 P1],B1',1);
lmiterm([1 6 1 0],-B*C1);

lmiterm([1 6 2 W],B1',1);
lmiterm([1 6 2 Cf11],-B,1);

lmiterm([1 6 6 0],-2*B*D1-C);
lmiterm([1 6 6 x],1,1);

lmiterm([1 7 1 R],tauM,A1);

lmiterm([1 7 6 R],tauM,B1);

lmiterm([1 7 7 R],-1,1);

lmiterm([1 8 1 0],M*C1);

lmiterm([1 8 2 Cf11],-M,1);

lmiterm([1 8 6 0],M*D1);

lmiterm([1 8 8 0],-I);

lmiterm([9 1 1 P1],1,A2,'s');
lmiterm([9 1 1 Q],1,1);
lmiterm([9 1 1 R],-1,1);

lmiterm([9 2 1 W],1,A2);
lmiterm([9 2 1 -Af22],1,1);

lmiterm([9 2 2 -Af22],1,1);
lmiterm([9 2 2 Af22],1,1);

lmiterm([9 3 1 -Bf22],L2',1);
lmiterm([9 3 1 R],1,1);
lmiterm([9 3 1 U],-1,1);

lmiterm([9 3 2 -Bf22],L2',1);

lmiterm([9 3 3 R],-2,1);
lmiterm([9 3 3 U],1,1);
lmiterm([9 3 3 -U],1,1);
lmiterm([9 3 3 O],g*L2',L2);

lmiterm([9 4 1 U],1,1);

lmiterm([9 4 3 R],1,1);
lmiterm([9 4 3 U],-1,1);

lmiterm([9 4 4 Q],-1,1);
lmiterm([9 4 4 R],-1,1);

lmiterm([9 5 1 -Bf22],-1,1);

lmiterm([9 5 2 -Bf22],1,1);

lmiterm([9 5 5 O],-1,1);

lmiterm([9 6 1 P1],B2',1);
lmiterm([9 6 1 0],-B*C2);

lmiterm([9 6 2 W],B2',1);
lmiterm([9 6 2 Cf22],-B,1);

lmiterm([9 6 6 0],-2*B*D2-C);
lmiterm([9 6 6 x],1,1);

lmiterm([9 7 1 R],0.1,A2);

lmiterm([9 7 6 R],0.1,B2);

lmiterm([9 7 7 R],-1,1);

lmiterm([9 8 1 0],M*C2);

lmiterm([9 8 2 Cf22],-M,1);

lmiterm([9 8 6 0],M*D2);

lmiterm([9 8 8 0],-I);

lmiterm([-3 1 1 P1],1,1);
lmiterm([-4 1 1 W],1,1);
lmiterm([-5 1 1 R],1,1);
lmiterm([-6 1 1 U],1,1);
lmiterm([-7 1 1 O],1,1);
lmiterm([-8 1 1 x],1,1);

lmisys=getlmis;
[tmin,xfeas]=feasp(lmisys);
P1=dec2mat(lmisys,xfeas,P1);
W=dec2mat(lmisys,xfeas,W);
R=dec2mat(lmisys,xfeas,R);
U=dec2mat(lmisys,xfeas,U);
O=dec2mat(lmisys,xfeas,O);
Af11=dec2mat(lmisys,xfeas,Af11);
Bf11=dec2mat(lmisys,xfeas,Bf11);
Cf11=dec2mat(lmisys,xfeas,Cf11);

Af22=dec2mat(lmisys,xfeas,Af22);
Bf22=dec2mat(lmisys,xfeas,Bf22);
Cf22=dec2mat(lmisys,xfeas,Cf22);

x=dec2mat(lmisys,xfeas,x);
% Df1=dec2mat(lmisys,xfeas,Df1);

Af1=Af11*W^-1;
Bf1=Bf11;
Cf1=Cf11*W^-1;
Df1=D1;

Af2=Af22*inv(W);
Bf2=Bf22;
Cf2=Cf22*W^-1;
Df2=D2;




 