% 矢量水听器DOA估计——单目标
clear ;close all;clc;
%% 信号
f1 = 60;
T = 400;
fs = 2048 ;
alpha = [-60] ;           % 方位角
theta = -180:1:180 ;
t = 0:1/fs:T-1/fs;
FL = 50;FH = 150 ;
noise = randn(1,T*fs);
b = fir1(128,[FL FH]/(fs/2)) ;
b2 = fir1(128,0.99) ;
noise2 = filter(b,1,noise) ;    % 带限噪声
xt =  0.1*sin(2*pi*f1*t) + 0.5*noise2 ;   % 构造声压通道信号：包含正弦信号与带限噪声
A = 3 ;
% 模拟矢量水听器的声压通道（p）和两个正交振速通道（x, y）
p = A*xt + filter(b2,1,randn(1,T*fs)) ;
x = A*xt .* cosd(alpha) + filter(b2,1,randn(1,T*fs)) ;
y = A*xt .* sind(alpha) + filter(b2,1,randn(1,T*fs)) ;
nbins = 90 ;    % 直方图这个地方取值不好设计，应该需要长时间积分

%% 处理
Nfft = 4096*2 ;
Len = 4096 ;
fLn = floor(FL/fs*Nfft)+1 ;
fHn = ceil(FH/fs*Nfft)+1 ;
% f1n = round(f1/fs*Nfft)+1 ;
% f2n = round(f2/fs*Nfft)+1 ;
Ktime = floor(T*fs/Len);
for k = 1 : Ktime
    % 当前数据段
    ptemp = p((k-1)*Len+1:k*Len) ;
    xtemp = x((k-1)*Len+1:k*Len) ;
    ytemp = y((k-1)*Len+1:k*Len) ;

    %% 平均声强器 —— 声压与振速的互谱实部
    Pf = fft(ptemp,Nfft)/Nfft ;
    Vx = fft(xtemp,Nfft)/Nfft ;
    Vy = fft(ytemp,Nfft)/Nfft ;

    Pvx = sum(Pf(fLn:fHn) .* conj(Vx(fLn:fHn))) ; % 共轭
    Pvy = sum(Pf(fLn:fHn) .* conj(Vy(fLn:fHn))) ;
    % 方位估计（取实部计算反正切）
    I1(k) = atand(real(Pvy)/real(Pvx));           % 各通道信号不存在相位差情况下取实部

    %% 直方图统计 —— 统计各频点角度估计的分布，取最大频次对应的角度
    Pvx2 = real(Pf(fLn:fHn) .* conj(Vx(fLn:fHn))) ;
    Pvy2 = real(Pf(fLn:fHn) .* conj(Vy(fLn:fHn))) ;
    Itemp = atand(Pvy2./Pvx2) ;
    % 统计角度直方图
    edges = linspace(-180, 180, nbins) ;
    N = histc(Itemp,edges); % 统计频次
    [~,pt] = max(N) ;
    I2_1(k) = (edges(pt)+edges(pt+1))/2;
    I2_2(k,:) = N;

    %% 加权直方图 —— 根据频谱幅度对角度估计进行加权统计，提高信噪比高的频点权重
    N2 = zeros(size(N)) ;
    Af = abs(Pf(fLn:fHn)) ;
    Afmean = mean(Af) ;
    AddN = round(Af./Afmean) ;  % 权值
    for k2 = 1 : length(Itemp) 
        for k3 = 1 : nbins
            if Itemp(k2)>edges(k3) && Itemp(k2)<edges(k3+1)
                N2(k3) = N2(k3) + AddN(k2);
                break ;
            end
        end
    end
    [~,v] = max(N2);
    I3_1(k) = (edges(v)+edges(v+1))/2;
    I3_2(k,:) = N2 ;

    %% DOA估计 —— CBF & MVDR & MUSIC
    Rx = [ptemp ;xtemp; ytemp]*[ptemp; xtemp; ytemp].'/Len ;    % 协方差矩阵
    invRx = inv(Rx+eye(3)*eps) ;    % 正则化逆矩阵
    % MUSIC算法：噪声子空间分解
    [V, D] = eig(Rx);
    [~, idx_sort] = sort(diag(D), 'descend');
    V = V(:, idx_sort);
    Rn = V(:, 2:3) * V(:, 2:3)'; % 噪声子空间投影矩阵
    for k3 = 1 : 361
        w = [1, cosd(theta(k3)), sind(theta(k3))]; % 导向矢量
        I4_cbf(k, k3) = w * Rx * w';               % CBF波束形成
        I4_mvdr(k, k3) = 1 / (w * invRx * w');     % MVDR波束形成
        I5_musuc(k, k3) = 1 / (w * Rn * w');       % MUSIC空间谱
    end
                 
end

%% 画图
figure ;
plot(I1,1:Ktime,'d','LineWidth',1.2);hold on;
plot(I2_1,1:Ktime,'p','LineWidth',1.2);hold on;
plot(I3_1,1:Ktime,'h','LineWidth',1.2) ;
legend('平均声强','直方图','加权直方图') ;
axis([-180 180 1 Ktime]);
figure;pcolor(linspace(-180,180,nbins),1:Ktime,I2_2);title('直方图法');shading interp;colormap(jet);
figure;pcolor(linspace(-180,180,nbins),1:Ktime,I3_2);title('加权直方图法');shading interp;colormap(jet);
figure;pcolor(theta,1:Ktime,10*log10(I4_cbf));title('CBF');shading interp;colormap(jet);
figure;pcolor(theta,1:Ktime,10*log10(I4_mvdr));title('MVDR');shading interp;colormap(jet);
figure;pcolor(theta,1:Ktime,10*log10(I5_musuc));title('MUSIC');shading interp;colormap(jet);

figure;
plot(theta,I4_cbf(1,:)/max(I4_cbf(1,:)));hold on
plot(theta,I4_mvdr(1,:)/max(I4_mvdr(1,:)));hold on
plot(theta,I5_musuc(1,:)/max(I5_musuc(1,:)));legend('CBF','MVDR','MUSIC') ;

figure;title('误差统计');
plot(I1-alpha);hold on ;
plot(I2_1-alpha);hold on ;
plot(I3_1-alpha);hold on ;
[~,v1] = max(I4_cbf.');
[~,v2] = max(I4_mvdr.');
[~,v3] = max(I5_musuc.');
plot(v1-alpha-1);hold on ;
plot(v2-alpha-1);hold on ;
plot(v3-alpha-1);hold on ;
legend('平均声强','直方图','加权直方图','CBF','MVDR','MUSIC');
disp(['平均声强  直方图  加权直方图  CBF  MVDR  MUSIC']);
Err=[mean(abs(I1-alpha)),mean(abs(I2_1-alpha)),mean(abs(I3_1-alpha)),mean(abs(v1-alpha-1)),mean(abs(v2-alpha-1)),mean(abs(v3-alpha-1))]


