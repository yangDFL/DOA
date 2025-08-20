%% Toeplitz重构
close all;clc; clear;
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
Theta=[20 30 115];                  %两个信号的入射角度
source_number = length(Theta);      %信元数
f = fs/N_x*(0:N_x-1);               %具体频点
fl = f0; fh = f0;

%% 生成信号
rng('default');
sig1 = cos(2*pi*f0* (t - d*(0:M-1)'*cosd(Theta(1)) / c));   % 目标信号
sig2 = cos(2*pi*f0* (t - d*(0:M-1)'*cosd(Theta(2)) / c));   
sig3 = cos(2*pi*f0* (t - d*(0:M-1)'*cosd(Theta(3)) / c)); 
noise = randn(M,length(sig1));                              % 均值为0，方差为1的高斯白噪声
snr = 0;
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
        %% TOP
        for n = 0:M-1
            temp = 0;
            for ii = 1:M-n
                temp = temp + R(ii,ii+n);
            end
            r(-n+M) = temp / (M-n);
            r(n+M) = conj(r(-n+M));
        end
        for ii = 1:M
            for jj = 1:M
                Rt1(ii,jj) = r(ii-jj+M);
            end
        end
        [V1,D1] = eig(Rt1);e1s=V1(:,M);       % 最大特征值对应的特征矢量
        Un1 = V1(:,1:M-source_number);

        %% MTOP
        for n = 0:M-1
            temp = 0;temp1=0;
            for ii = 1:M-n
                temp = temp + abs(R(ii,ii+n));
                temp1 = temp1 + R(ii,ii+n);
            end
            phase = atan2(imag(temp1),real(temp1));
            r1(-n+M) = temp / (M-n) * exp(1j*phase);
            r1(n+M) = conj(r1(-n+M));
        end
        for ii = 1:M
            for jj = 1:M
                Rt2(ii,jj) = r1(ii-jj+M);
            end
        end
        [V2,D2] = eig(Rt2);es2=V2(:,M);       % 最大特征值对应的特征矢量
        Un2 = V2(:,1:M-source_number);

        p_cbf(i,j) = atheta'*x(:,j);
        p_top(i,j) = 1 / (atheta'*Un1*Un1'*atheta);
        p_mtop(i,j) = 1 / (atheta'*Un2*Un2'*atheta);

    end
end

P_cbf = sum(p_cbf.' .* conj(p_cbf.')) / length(find(f>=fl & f<=fh)); 
P_top = sum(p_top.' .* conj(p_top.')) / length(find(f>=fl & f<=fh)); 
P_mtop = sum(p_mtop.' .* conj(p_mtop.')) / length(find(f>=fl & f<=fh)); 

%% 画图
figure;
plot(theta,10*log10(P_cbf/max(P_cbf)));
hold on;
plot(theta,10*log10(P_top/max(P_top)));
plot(theta,10*log10(P_mtop/max(P_mtop)));
xlabel('角度/°');ylabel('幅度/dB');



