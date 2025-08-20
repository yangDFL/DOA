%-----------------------------------------%
%           xulingji 2011.7.7             %
% 双边相关变换（TCT）估计宽带相干信号源方位 %
%-----------------------------------------%
clear all
close all
clc
tic

jay=sqrt(-1);
M=16;               %阵元数
N=15000;            %仿真产生的快拍数
N_use=12800;        %实际用的快拍数
n1=1000;            %数据起始点          
f0=2000;            %以f0为聚焦频率
f_l=3*f0/5;         %宽带信号上边界频率
f_u=f0;             %宽带信号下边界频率
delta_F=40;
F=f_l:delta_F:f_u;
B=(f_u-f_l);        %绝对带宽
B_delta=(f_u-f_l)/((f_l+f_u)/2);  %相对带宽（一般取50%）%%%这样子难道频率上下界都要有规律？
Fs=5*f0;
c=1500;
lambda=c/f0;        %f0对应的波长
d=1/2*lambda;       %阵元间距
rsj=[-1.5 1];       
Num=length(rsj);
Q1=-90;
Q2=90;
deta_theta=0.02;
Angles=Q1:deta_theta:Q2;
z=0:M-1;
SNR=15; 

L=100;               
k1=-L:L;
index=1;
nfft=256;
delta_f=Fs/nfft;
f1=round(f_l/delta_f+1);    %频率对应于FFT中的点数
f2=round(f_u/delta_f+1);
freq=f1:f2;
freq_hz=(freq-1)*Fs/nfft;   %实际分析的频率

Q=512;                        
wn1=800*2/Fs;wn2=2200*2/Fs;   
b=fir1(Q,[wn1,wn2]);
no1=randn(1,N);
no2=randn(1,N);
no3=randn(M,N);
s1=filter(b,1,no1);
s2=filter(b,1,no1);
noise=filter(b,1,no3,[],2);  %噪声
t0=10*1/Fs;
for i=1:M
    tao1(i)=d*sin(rsj(1)*pi/180)/c*z(i)*Fs;        %入射信号1延迟点数 
    tao2(i)=d*sin(rsj(2)*pi/180)/c*z(i)*Fs+t0*Fs;  %入射信号2延迟点数 
end
for j=1:M
    Coef1=sinc(k1+tao1(j));
    x1(j,:)=filter(Coef1,1,s1);
    Coef2=sinc(k1+tao2(j));
    x2(j,:)=filter(Coef2,1,s2);
end
S=x1+x2+(10^(-SNR/20))*noise;
S=S(:,n1:n1+N_use-1);

% % % % % % % % % % % % % % % % % % % % % % % % % % % 
Section=1;                      %分段数的起始值
Rk=zeros(M,M,length(freq));
while index+nfft-1 < size(S,2)  %构造出频域相关矩阵
    X=fft(S(:,index:index+nfft-1),nfft,2);
%     Nf=fft(No(:,index:index+nfft-1),nfft,2);
    for ii=1:length(freq)
        Rk(:,:,ii)=Rk(:,:,ii)+X(:,freq(ii))*X(:,freq(ii))';
%         Rn(:,:,ii)=Rn(:,:,ii)+Nf(:,freq(ii))*Nf(:,freq(ii))';
    end
    index=index+nfft;
    Section=Section+1;
end
rsj_predict=[ -1.5 -1.3 0.99 1.06];
mu=zeros(M,1);
 sj=zeros(length(rsj_predict));
for jj=1:length(freq)      %寻找出聚焦矩阵
      A=exp(jay*2*pi*freq_hz(ii)/c*z'*d*sin(rsj_predict*pi/180));
%       R=1/Section*reshape(Rk(:,:,jj),M,M);
      [Ur Dr Vr]=svd(Rk(:,:,jj));
      Vj(:,:,jj)=Vr(:,1:Num);
      sigma=mean(diag(Dr(Num+1:end,Num+1:end)));%%对应的噪声
      Pj=Rk(:,:,jj)-sigma*eye(M);
      [Upj Dpj Vpj]=svd(Pj);
      mu=mu+diag(Dpj);
      sj=sj+inv(A'*A)*A'*Pj*A*inv(A'*A);
end
  s0=1/length(freq)*sj;
  mu=1/length(freq)*mu;
 for k=1:length(F)            %选择参考频率fc
    A0=exp(jay*2*pi*F(k)/c*z'*d*sin(rsj_predict*pi/180));
    P0=A0*s0*A0';
    [Up Dp Vp]=svd(P0);
    Dp0=diag(Dp);
    SUM=0;
    for kk=1:Num
        SUM=SUM+(Dp0(kk)-mu(kk)).^2;
    end
    Result(k)=SUM;
 end
[Y I]=min(Result);
fc=F(1)+(I-1)*delta_F;
A0=exp(jay*2*pi*fc/c*z'*d*sin(rsj_predict*pi/180));
P0=A0*s0*A0';
[U0 D0 V0]=svd(P0);
Ry=zeros(M);
for jj=1:length(freq)    %变换信号相关矩阵
    T=V0(:,1:Num)*Vj(:,:,jj)';
    Ry=Ry+T*Rk(:,:,jj)*T';
end
Ry=1/length(freq)*Ry;
a=exp(jay*2*pi*fc/c*z'*d*sin(Angles*pi/180));
[V D]=svd(Ry);        %广义特征分解后常规MUSIC估计方位角
En=V(:,Num+1:end);
P_music=diag(abs((a'*a)./(a'*En*En'*a)));
P_CSM=10*log10(P_music/max(P_music));
figure,
plot(Angles,P_CSM,'r','linewidth',1)
set(gca,'fontsize',14)
xlabel('方位 /度')
ylabel('方位谱 /dB')
title('TCT-CSM法估计宽带信号方位')
grid on
box on
%如果有必要迭代初始方位角，提高工作性能
toc



