clear;
clc;
%%
Sensor=8;                                     %天线阵元数
bw1=1e7;                                      % 信号带宽
bw2=1e6;
bw3=51e6;
f1=1e8;
f2=1e7;
f3=1e7;
T1=1e-2;                                      %信号脉冲宽度  
L=1024;                                       %采样点数
dt1=T1/L;                                     %采样间隔
gc=3e8;                                       %光速
snr=-500;                                       %信躁比
degrad=pi/180;
p=3;                                          %信号源个数 
angle=1;           
angle2=2;                                    %信号源角度
angle3=-4;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%构造采样序列
t1=(0:L-1)*dt1;                               %时域采样点
sig1=exp(j*pi*bw1*t1.*t1/T1+j*2*pi*f1*t1);    %产生宽带信号
sig2=exp(j*pi*bw2*t1.*t1/T1+j*2*pi*f2*t1);  
sig3=exp(j*pi*bw3*t1.*t1/T1+j*2*pi*f3*t1);  
sig1=awgn(sig1,snr);
sig2=awgn(sig2,snr);                          %加白噪音
sig3=awgn(sig3,snr);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%宽带信号频率
F=f1+bw1/(T1)*t1;                           
f1=F(128);
f2=F(2*128);
f3=F(3*128);
f4=F(4*128);
f5=F(5*128);
f6=F(6*128);
f7=F(7*128);
f8=F(8*128);
F=[f1,f2,f3,f4,f5,f6,f7,f8];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%构造采样数据矩阵
Lambda=gc/f8;                                 %波长
d=0.5*Lambda;                                 %阵元间距离
tmp2=d*[0:Sensor-1 ]';            

v11=sig1(1:128);
v12=sig1(129:256);
v13=sig1(257:384);
v14=sig1(385:512);
v15=sig1(513:640);
v16=sig1(641:768);
v17=sig1(769:896);
v18=sig1(897:1024);
v21=sig2(1:128);
v22=sig2(129:256);
v23=sig2(257:384);
v24=sig2(385:512);
v25=sig2(513:640);
v26=sig2(641:768);
v27=sig2(769:896);
v28=sig2(897:1024);
v31=sig3(1:128);
v32=sig3(129:256);
v33=sig3(257:384);
v34=sig3(385:512);
v35=sig3(513:640);
v36=sig3(641:768);
v37=sig3(769:896);
v38=sig3(897:1024);

s1=[(v11);(v21);v31];
s2=[(v12);(v22);v32];
s3=[(v13);(v23);v33];
s4=[(v14);(v24);v34];
s5=[(v15);(v25);v35];
s6=[(v16);(v26);v36];
s7=[(v17);(v27);v37];
s8=[(v18);(v28);v38];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%构造方向矩阵
A1=[exp(-j*2*pi*f1/gc*tmp2*sin(angle*degrad)),exp(-j*2*pi*f1/gc*tmp2*sin(angle2*degrad)),exp(-j*2*pi*f1/gc*tmp2*sin(angle3*degrad))]; 
A2=[exp(-j*2*pi*f2/gc*tmp2*sin(angle*degrad)),exp(-j*2*pi*f2/gc*tmp2*sin(angle2*degrad)),exp(-j*2*pi*f2/gc*tmp2*sin(angle3*degrad))]; 
A3=[exp(-j*2*pi*f3/gc*tmp2*sin(angle*degrad)),exp(-j*2*pi*f3/gc*tmp2*sin(angle2*degrad)),exp(-j*2*pi*f3/gc*tmp2*sin(angle3*degrad))]; 
A4=[exp(-j*2*pi*f4/gc*tmp2*sin(angle*degrad)),exp(-j*2*pi*f4/gc*tmp2*sin(angle2*degrad)),exp(-j*2*pi*f4/gc*tmp2*sin(angle3*degrad))]; 
A5=[exp(-j*2*pi*f5/gc*tmp2*sin(angle*degrad)),exp(-j*2*pi*f5/gc*tmp2*sin(angle2*degrad)),exp(-j*2*pi*f5/gc*tmp2*sin(angle3*degrad))]; 
A6=[exp(-j*2*pi*f6/gc*tmp2*sin(angle*degrad)),exp(-j*2*pi*f6/gc*tmp2*sin(angle2*degrad)),exp(-j*2*pi*f6/gc*tmp2*sin(angle3*degrad))]; 
A7=[exp(-j*2*pi*f7/gc*tmp2*sin(angle*degrad)),exp(-j*2*pi*f7/gc*tmp2*sin(angle2*degrad)),exp(-j*2*pi*f7/gc*tmp2*sin(angle3*degrad))]; 
A8=[exp(-j*2*pi*f8/gc*tmp2*sin(angle*degrad)),exp(-j*2*pi*f8/gc*tmp2*sin(angle2*degrad)),exp(-j*2*pi*f8/gc*tmp2*sin(angle3*degrad))]; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%接受信号
SNR=10;
x1=A1*s1;
x2=A2*s2;
x3=A3*s3;
x4=A4*s4;
x5=A5*s5;
x6=A6*s6;
x7=A7*s7;
x8=A8*s8;
x1=awgn(x1,SNR);
x2=awgn(x2,SNR);
x3=awgn(x3,SNR);
x4=awgn(x4,SNR);
x5=awgn(x5,SNR);                              %加白噪音
x6=awgn(x6,SNR);
x7=awgn(x7,SNR);
x8=awgn(x8,SNR);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FFT
X1=fft(x1.');
X2=fft(x2.');
X3=fft(x3.');
X4=fft(x4.');
X5=fft(x5.'); 
X6=fft(x6.');
X7=fft(x7.');
X8=fft(x8.');
X1=X1.';
X2=X2.';
X3=X3.';
X4=X4.';
X5=X5.';
X6=X6.';
X7=X7.';
X8=X8.';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%计算协方差矩阵
R1=X1*X1';
R2=X2*X2';
R3=X3*X3';
R4=X4*X4';
R5=X5*X5';
R6=X6*X6';
R7=X7*X7';
R8=X8*X8';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%计算去噪后的协方差矩阵 
temp=eig(R1);                                    
temp=sort(temp);
P1=R1-[sum(temp)-temp(8)-temp(7)-temp(6)]/(8-p)*eye(8,8);
temp=eig(R2);
temp=sort(temp);
P2=R2-[sum(temp)-temp(8)-temp(7)-temp(6)]/(8-p)*eye(8,8);
temp=eig(R3);
temp=sort(temp);
P3=R3-[sum(temp)-temp(8)-temp(7)-temp(6)]/(8-p)*eye(8,8);
temp=eig(R4);
temp=sort(temp);
P4=R4-[sum(temp)-temp(8)-temp(7)-temp(6)]/(8-p)*eye(8,8);
temp=eig(R5);
temp=sort(temp);
P5=R5-[sum(temp)-temp(8)-temp(7)-temp(6)]/(8-p)*eye(8,8);
temp=eig(R6);
temp=sort(temp);
P6=R6-[sum(temp)-temp(8)-temp(7)-temp(6)]/(8-p)*eye(8,8);
temp=eig(R7);
temp=sort(temp);
P7=R7-[sum(temp)-temp(8)-temp(7)-temp(6)]/(8-p)*eye(8,8);
temp=eig(R8);
temp=sort(temp);
P8=R8-[sum(temp)-temp(8)-temp(7)-temp(6)]/(8-p)*eye(8,8);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%奇异值分解
[U1,S1,V1]=svd(P1);
[U2,S2,V2]=svd(P2);
[U3,S3,V3]=svd(P3);
[U4,S4,V4]=svd(P4);
[U5,S5,V5]=svd(P5);
[U6,S6,V6]=svd(P6);
[U7,S7,V7]=svd(P7);
[U8,S8,V8]=svd(P8);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%选择中心频率          
S=[sum(S1);sum(S2);sum(S3);sum(S4);sum(S5);sum(S6);sum(S7);sum(S8)];
S_=sum(S);
sum=zeros(1,8);
for j=1:8
    for i=1:8
        delta=S(i,j)-S_(1,i)/8;                                    
        sum(1,j)=sum(1,j)+delta^2;
    end
end
[h,j]=min(sum);   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 构造各频率点聚焦矩阵
[Q1,D1]=eig(P1);
[Q2,D2]=eig(P2);
[Q3,D3]=eig(P3);
[Q4,D4]=eig(P4);
[Q5,D5]=eig(P5);
[Q6,D6]=eig(P6);
[Q7,D7]=eig(P7);
[Q8,D8]=eig(P8);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%选择聚焦矩阵
Q0=[Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8];
Q_=Q0(1+(j-1)*64:j*64);
Q=reshape(Q_,8,8);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%聚焦变换
T1=Q*Q1';      
T2=Q*Q2';
T3=Q*Q3';
T4=Q*Q4';
T5=Q*Q5';
T6=Q*Q6';
T7=Q*Q7';
T8=Q*Q8';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MUSIC算法估计信号入射方向
RY=(T1*X1)*(T1*X1)'+(T2*X2)*(T2*X2)'+(T3*X3)*(T3*X3)'+(T4*X4)*(T4*X4)'+(T5*X5)*(T5*X5)'+(T6*X6)*(T6*X6)'+(T7*X7)*(T7*X7)'+(T8*X8)*(T8*X8)';
[U,S11,V]=svd(RY);
Vs=U(:,1:p);
Vu=U(:,p+1:Sensor);
i= sqrt(-1);
th2=[-90:0.1:90]';
T=length(th2);
tmp=-i* 2*pi*F(j)/gc *sin(th2'* degrad);
a2=tmp2*tmp;
a2=exp(a2);
num=diag(a2'*a2);                             %分子
Ena=Vu'* a2;
den=diag(Ena'* Ena);                          %分母
doa=(num./den);                               %空间谱函数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%作图
for i=1:T
    doa1(i)=20*log10(abs(doa(i))/max(abs(doa)));
end
plot(th2, doa1,'-');
title('TCT算法');
xlabel('角度/°');                               %横轴
ylabel('空间谱/dB');                             %纵轴
% axis([-90 90 0.1 1e6]);
grid on;
hold on;

