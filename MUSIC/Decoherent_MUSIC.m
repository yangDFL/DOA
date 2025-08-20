%% 解相干MUSIC——空间平滑、矩阵重构
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
Theta=[20 45 115];                  %两个信号的入射角度
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
N = 14;                         % 平滑后子阵数量
sub_M = M-N+1;                  % 子阵阵元数

for i = 1 : length(theta)
    for j = find(f>=fl & f<=fh)
        R = x(:,j) * x(:,j)' / N_x;
        atheta = exp(-1j*(0:M-1)'*2*pi*d*f(j)*cosd(theta(i))/c);
        atheta1 = exp(-1j*(0:N-1)'*2*pi*d*f(j)*cosd(theta(i))/c);
        %%% 前向平滑
        Rx1 = zeros(N,N);
        for ii = 1:sub_M   %第1个子阵列到第N个子阵列
            Rx1 = Rx1 + R(ii:ii+N-1, ii:ii+N-1);%向前滑动，注意每个子阵列是k个阵元，不是k-1个
        end
        R1 = Rx1 / sub_M;
        [V1,D1] = eig(R1);
        Un1 = V1(:,1:end-source_number);
        %%% 后向平滑
        J = fliplr(eye(M));
        RR = J*conj(R)*J;
        Rx3 = zeros(N,N);
        for ii = 1:sub_M
            Rx3 = Rx3 + RR(ii:ii+N-1, ii:ii+N-1);
        end
        R4 = Rx3 / sub_M;
        [V4,D4] = eig(R4);
        Un4 = V4(:,1:end-source_number);
        %%% 双向平滑
        J = fliplr(eye(M));                 % 翻转对角阵
        Rx = (R + J*conj(R)*J)/2;           % 前后向平滑协方差矩阵
        Rx2 = zeros(N,N);
        for ii = 1:sub_M                    % 第1个子阵列到第N个子阵列
            Rx2 = Rx2 + Rx(ii:ii+N-1,ii:ii+N-1);%向前滑动，注意每个子阵列是k个阵元，不是k-1个
        end
        R2 = Rx2 / sub_M;
        [V2,D2] = eig(R2);
        Un2 = V2(:,1:end-source_number);
        %%% 传统MUSIC
        [V,D] = eig(R);es=V(:,M);       % 最大特征值对应的特征矢量
        Un = V(:,1:M-source_number);
        %%% 重构矩阵
        for ii=1:m
            Y(ii,:)=es(ii:ii+p-1);
        end
        [V3,D3]=svd(Y); En = V3(:,source_number+1:m);
        b_theta = exp(-1j*(0:m-1)'*2*pi*d*f(j)*cosd(theta(i))/c);
        
        p_re(i,j) = 1 / (b_theta'*En*En'*b_theta);
        p_fmusic(i,j) = 1 / (atheta1'*Un1*Un1'*atheta1);
        p_bmusic(i,j) = 1 / (atheta1'*Un4*Un4'*atheta1);
        p_fbmusic(i,j) = 1 / (atheta1'*Un2*Un2'*atheta1);
        p_music(i,j) = 1 / (atheta'*Un*Un'*atheta);
        p_cbf(i,j) = atheta'*x(:,j);
        
    end
end
P_cbf = sum(p_cbf.' .* conj(p_cbf.')) / length(find(f>=fl & f<=fh)); 
P_music = sum(p_music.' .* conj(p_music.')) / length(find(f>=fl & f<=fh));  
P_fmusic = sum(p_fmusic.' .* conj(p_fmusic.')) / length(find(f>=fl & f<=fh));  
P_bmusic = sum(p_bmusic.' .* conj(p_bmusic.')) / length(find(f>=fl & f<=fh));
P_fbmusic = sum(p_fbmusic.' .* conj(p_fbmusic.')) / length(find(f>=fl & f<=fh));
P_re = sum(p_re.' .* conj(p_re.')) / length(find(f>=fl & f<=fh));

%% 画图
figure;
plot(theta,10*log10(P_cbf/max(P_cbf)));
hold on;
plot(theta,10*log10(P_music/max(P_music)));
plot(theta,10*log10(P_fmusic/max(P_fmusic)));
plot(theta,10*log10(P_bmusic/max(P_bmusic)));
plot(theta,10*log10(P_fbmusic/max(P_fbmusic)));
plot(theta,10*log10(P_re/max(P_re)));
xlabel('角度/°');ylabel('幅度/dB');
legend('CBF','MUSIC','Smooth-F\_MUSIC','Smooth-B\_MUSIC','Smooth-FB\_MUSIC','Vector\_Refactoring');





