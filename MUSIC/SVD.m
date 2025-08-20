%% 矢量奇异值（SVD）
close all; clc; clear;
%%
c = 1500;
M = 20;                             %阵元数
m = 10;
p = M-m+1;
f0 = 500;
N_x = 4096;                         %信号长度
T = 1;
t = 0:1/N_x:T-1/N_x;                %信号时间
fs = N_x;                           %快拍数
l = c / f0;                         %信号波长  
d = 0.5*l;                          %阵元间距
Theta=[20 25 100];                  %两个信号的入射角度
source_number = length(Theta);      %信元数
f = fs/N_x*(0:N_x-1);               %具体频点
fl = f0; fh = f0;
%% 生成信号
rng('default');
sig1 = cos(2*pi*f0* (t - d*(0:M-1)'*cosd(Theta(1)) / c));   % 目标信号
sig2 = cos(2*pi*f0* (t - d*(0:M-1)'*cosd(Theta(2)) / c));   
sig3 = cos(2*pi*f0* (t - d*(0:M-1)'*cosd(Theta(3)) / c)); 
noise = randn(M,length(sig1));                              % 均值为0，方差为1的高斯白噪声
snr = 10;
SNR1 = 10^(snr / 10);  SNR2 = 10^(snr / 10); SNR3 = 10^(snr/10);
P1 = mean(sig1.^2);  P2 = mean(sig2.^2); P3 = mean(sig3.^2);                 % 计算信号功率
PN1 = P1 / SNR1;       PN2 = P2 / SNR2;  PN3 = P3 / SNR3;                  % 计算所需的噪声功率
noise_scale1 = sqrt(PN1 / mean(noise.^2));                  % 计算噪声缩放因子
noise_scale2 = sqrt(PN2 / mean(noise.^2));
noise_scale3 = sqrt(PN3 / mean(noise.^2));
sig1 = sig1 + noise_scale1*noise;
sig2 = sig2 + noise_scale2*noise;
sig3 = sig3 + noise_scale3*noise;
Sig = sig1 + sig2 + sig3;
%% DOA估计
theta = linspace(0,180,180);
x = fft(Sig').';
for i = 1 : length(theta)
    for j = find(f>=fl & f<=fh)
        R = x(:,j) * x(:,j)' / N_x;
        atheta = exp(-1j*(0:M-1)'*2*pi*d*f(j)*cosd(theta(i))/c);
        atheta1 = exp(-1j*(0:m-1)'*2*pi*d*f(j)*cosd(theta(i))/c);
        [U,S,V] = svd(R);
        %% 方式1 ESVD
        e1 = U(:,1);    % 最大特征值对应的特征矢量
        %% 方式2 DSVD
        X1 = x(:,j)*x(1,j)' / N_x;
        e2 = X1;
        %% 方式3 PSVD
        % 需要先确定信号的大致方向
        theta0 = 20;
        atheta0 = exp(-1j*(0:M-1)'*2*pi*d*f(j)*cosd(theta0)/c);
        Bk = atheta0'*x(:,j) / M;
        X2 = x(:,j)*Bk / N_x;
        e3 = X2;

        for ii = 1:m
            Y1(ii,:) = e1(ii:ii+p-1);
            Y2(ii,:) = e2(ii:ii+p-1);
            Y3(ii,:) = e3(ii:ii+p-1);
        end
        [U1,S1,V1] = svd(Y1);   Un1 = U1(:,source_number+1:end);
        [U2,V2,S2] = svd(Y2);   Un2 = U2(:,source_number+1:end);
        [U3,V3,S3] = svd(Y3);   Un3 = U3(:,source_number+1:end);

        p_cbf(i,j) = atheta'*x(:,j);
        p1(i,j) = 1 / (atheta1'*Un1*Un1'*atheta1);
        p2(i,j) = 1 / (atheta1'*Un2*Un2'*atheta1);
        p3(i,j) = 1 / (atheta1'*Un3*Un3'*atheta1);
    end
end
P_cbf = sum(p_cbf.' .* conj(p_cbf.')) / length(find(f>=fl & f<=fh)); 
P1 = sum(p1.' .* conj(p1.')) / length(find(f>=fl & f<=fh)); 
P2 = sum(p2.' .* conj(p2.')) / length(find(f>=fl & f<=fh)); 
P3 = sum(p3.' .* conj(p3.')) / length(find(f>=fl & f<=fh)); 

%% 画图
figure;
plot(theta,10*log10(P_cbf/max(P_cbf)));
hold on;
plot(theta,10*log10(P1/max(P1)));
plot(theta,10*log10(P2/max(P2)));
plot(theta,10*log10(P3/max(P3)));
xlabel('角度/°');ylabel('幅度/dB');



