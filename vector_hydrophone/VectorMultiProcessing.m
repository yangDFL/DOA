%% 矢量水听器多方位角DOA估计程序
clear all; close all; clc;

%% 参数设置
fs = 2048;                 % 采样频率
T = 400;                   % 信号时长
f1 = 60;                   % 信号频率
f2 = 120;
alpha = [-40, 80];         % 目标方位角（两个方向）
SNR = -10;                  % 信噪比
Nfft = fs;               % FFT点数
Len = Nfft;                % 每次处理的数据长度
Ktime = floor(T*fs/Len);   % 处理次数
nbins = 360;              % 直方图区间数
theta = linspace(-180,180,nbins+1);            % 方位角搜索范围

%% 滤波器设计
% 带通滤波器(50-150Hz)
[b, a] = fir1(64, [50 150]/(fs/2));
% 低通滤波器(截止频率接近奈奎斯特频率)
[b2, a2] = fir1(64, 0.99);

%% 频率范围设置
fL = 50;                  % 低频截止频率
fH = 150;                 % 高频截止频率
fLn = round(fL*Nfft/fs);  % 低频对应FFT索引
fHn = round(fH*Nfft/fs);  % 高频对应FFT索引
f_axis = (0:Nfft-1)*fs/Nfft; % 频率轴

%% 初始化结果存储
theta_Iavg = zeros(Ktime, 1);      % 平均声强法估计角度
theta_hist = zeros(Ktime, 2);      % 直方图法估计角度(两个方位角)
theta_weight = zeros(Ktime, 2);    % 加权直方图法估计角度(两个方位角)

%% 主循环 - 分时段处理数据
for k = 1:Ktime
    % 显示进度
    if mod(k, 50) == 0
        fprintf('Processing segment %d of %d\n', k, Ktime);
    end
    
    % 生成时间序列
    t = ((k-1)*Len+1 : k*Len)/fs;
    
    %% 生成信号
    % 信号由两个不同方位角的分量组成
    xt1 = sin(2*pi*f1*t);  % 第一个信号分量
    xt2 = sin(2*pi*f2*t);  % 第二个信号分量
    
    % 带限噪声
    nt = randn(size(t));
    nt = filter(b, a, nt);
    nt = nt / std(nt);
    
    % 合成信号
    xt = xt1 + xt2 + 10^(-SNR/20) * nt;
    
    % 生成矢量水听器信号
    p = 3*xt + randn(size(xt)) * 10^(-SNR/20);  % 声压通道
    x = 3*xt1 * cosd(alpha(1)) + 3*xt2 * cosd(alpha(2)) + randn(size(xt)) * 10^(-SNR/20);  % x方向振速
    y = 3*xt1 * sind(alpha(1)) + 3*xt2 * sind(alpha(2)) + randn(size(xt)) * 10^(-SNR/20);  % y方向振速
    
    % 滤波
    p = filter(b, a, p);
    x = filter(b, a, x);
    y = filter(b, a, y);
    
    % 零填充到Nfft长度
    p_fft = fft(p, Nfft) / Nfft;
    x_fft = fft(x, Nfft) / Nfft;
    y_fft = fft(y, Nfft) / Nfft;
    
    % 计算声压谱
    Ap = abs(p_fft);
    
    % 计算互谱
    Pvx = sum(p_fft(fLn:fHn) .* conj(x_fft(fLn:fHn)));
    Pvy = sum(p_fft(fLn:fHn) .* conj(y_fft(fLn:fHn)));
    
    %% 平均声强法
    theta_Iavg(k) = atand(real(Pvy)/real(Pvx));
    % theta_Iavg(k) = mod(Itemp + 180, 360);
        
    %% 直方图法
    Pvx2 = real(p_fft(fLn:fHn) .* conj(x_fft(fLn:fHn)));
    Pvy2 = real(p_fft(fLn:fHn) .* conj(y_fft(fLn:fHn)));
    Itemp1 = atand(Pvy2./Pvx2);
    
    % 角度归一化到-180~180度
    for i = 1:length(Itemp1)
        if Itemp1(i) < -180
            Itemp1(i) = Itemp1(i) + 360;
        elseif Itemp1(i) > 180
            Itemp1(i) = Itemp1(i) - 360;
        end
    end
    
    % 统计直方图
    count = histcounts(Itemp1, linspace(-180, 180, nbins+1));
    [~, idx] = sort(count, 'descend');
    
    % 找出两个峰值对应的角度区间
    bin_edges = linspace(-180, 180, nbins+1);
    bin_centers = (bin_edges(1:end-1) + bin_edges(2:end))/2;
    
    % 选择前两个峰值对应的角度
    theta_hist(k, 1) = bin_centers(idx(1));
    theta_hist(k, 2) = bin_centers(idx(2));
    I1(k,:) = count;
    
    %% 加权直方图法
    Af = Ap(fLn:fHn);
    AddN = Af / mean(Af);
    AddN = AddN / max(AddN);
    
    % 加权统计
    count_weight = zeros(1, nbins);
    for i = 1:length(Itemp1)
        bin_idx = find(Itemp1(i) <= bin_edges, 1) - 1;
        if bin_idx > 0 && bin_idx <= nbins
            count_weight(bin_idx) = count_weight(bin_idx) + AddN(i);
        end
    end
    
    [~, idx_weight] = sort(count_weight, 'descend');
    theta_weight(k, 1) = bin_centers(idx_weight(1));
    theta_weight(k, 2) = bin_centers(idx_weight(2));
    I2(k,:) = count_weight;

    %% DOA估计 —— CBF & MVDR & MUSIC
    Rx = [p ;x; y]*[p; x; y].'/Len ;    % 协方差矩阵
    invRx = inv(Rx+eye(3)*eps) ;    % 正则化逆矩阵
    % MUSIC算法：噪声子空间分解
    [V, D] = eig(Rx);
    Rn = V(:, 1) * V(:, 1)'; % 噪声子空间投影矩阵
    for k3 = 1 : length(theta)
        w = [1, cosd(theta(k3)), sind(theta(k3))]; % 导向矢量
        I4_cbf(k, k3) = w * Rx * w';               % CBF波束形成
        I4_mvdr(k, k3) = 1 / (w * invRx * w');     % MVDR波束形成
        I5_musuc(k, k3) = 1 / (w * Rn * w');       % MUSIC空间谱
    end
end

%% 结果可视化
% 时间轴
time_axis = (1:Ktime) * Len / fs;

% 平均声强法结果
figure;
plot(time_axis, theta_Iavg', 'b-', 'LineWidth', 1.5);
hold on;
plot(time_axis, alpha(1)*ones(size(time_axis)), 'r--', 'LineWidth', 1.5);
plot(time_axis, alpha(2)*ones(size(time_axis)), 'r--', 'LineWidth', 1.5);
grid on;
title('平均声强法 DOA 估计');
xlabel('时间 (s)');
ylabel('角度 (度)');
legend('估计角度', '真实角度1', '真实角度2');
ylim([-180, 180]);

% 直方图法结果
figure;
plot(time_axis, theta_hist(:,1), 'b-', 'LineWidth', 1.5);
hold on;
plot(time_axis, theta_hist(:,2), 'g-', 'LineWidth', 1.5);
plot(time_axis, alpha(1)*ones(size(time_axis)), 'r--', 'LineWidth', 1.5);
plot(time_axis, alpha(2)*ones(size(time_axis)), 'r--', 'LineWidth', 1.5);
grid on;
title('直方图法 DOA 估计');
xlabel('时间 (s)');
ylabel('角度 (度)');
legend('估计角度1', '估计角度2', '真实角度1', '真实角度2');
ylim([-180, 180]);

% 加权直方图法结果
figure;
plot(time_axis, theta_weight(:,1), 'b-', 'LineWidth', 1.5);
hold on;
plot(time_axis, theta_weight(:,2), 'g-', 'LineWidth', 1.5);
plot(time_axis, alpha(1)*ones(size(time_axis)), 'r--', 'LineWidth', 1.5);
plot(time_axis, alpha(2)*ones(size(time_axis)), 'r--', 'LineWidth', 1.5);
grid on;
title('加权直方图法 DOA 估计');
xlabel('时间 (s)');
ylabel('角度 (度)');
legend('估计角度1', '估计角度2', '真实角度1', '真实角度2');
ylim([-180, 180]);

figure;pcolor(theta,1:Ktime,10*log10(I4_cbf));title('CBF');shading interp;colormap(jet);
figure;pcolor(theta,1:Ktime,10*log10(I4_mvdr));title('MVDR');shading interp;colormap(jet);
figure;pcolor(theta,1:Ktime,10*log10(I5_musuc));title('MUSIC');shading interp;colormap(jet);


%% 方位历程图
figure;pcolor(linspace(-180,180,nbins),1:Ktime,I1);
title('直方图法');shading interp;colormap(jet);
figure;pcolor(linspace(-180,180,nbins),1:Ktime,I2);
title('加权直方图法');shading interp;colormap(jet);

%% 误差分析
% 计算每个估计值与最近真实值的误差
err_Iavg = zeros(Ktime, 1);
err_hist = zeros(Ktime, 2);
err_weight = zeros(Ktime, 2);
err_cbf = zeros(Ktime, 2);
err_mvdr = zeros(Ktime, 2);
err_music = zeros(Ktime, 2);

for k = 1:Ktime
    % 平均声强法误差(选择最近的真实角度)
    err1 = abs(theta_Iavg(k) - alpha(1));
    err2 = abs(theta_Iavg(k) - alpha(2));
    err_Iavg(k) = min(err1, err2);

    % 直方图法误差(匹配每个估计角度到最近的真实角度)
    for i = 1:2
        err1 = abs(theta_hist(k,i) - alpha(1));
        err2 = abs(theta_hist(k,i) - alpha(2));
        err_hist(k,i) = min(err1, err2);
    end

    % 加权直方图法误差
    for i = 1:2
        err1 = abs(theta_weight(k,i) - alpha(1));
        err2 = abs(theta_weight(k,i) - alpha(2));
        err_weight(k,i) = min(err1, err2);
    end

    % CBF & MVDR & MUSIC 误差
    [~,v1] = max(I4_cbf.');
    [~,v2] = max(I4_mvdr.');
    [~,v3] = max(I5_musuc.');
    for i = 1:2
        err1 = abs(theta(v1(1)) - alpha(1));
        err2 = abs(theta(v1(2)) - alpha(2));
        err_cbf(k,i) = min(err1, err2);

        err1 = abs(theta(v2(1)) - alpha(1));
        err2 = abs(theta(v2(2)) - alpha(2));
        err_mvdr(k,i) = min(err1, err2);

        err1 = abs(theta(v3(1)) - alpha(1));
        err2 = abs(theta(v3(2)) - alpha(2));
        err_music(k,i) = min(err1, err2);
    end
    
end

% 计算平均误差
avg_err_Iavg = mean(err_Iavg);
avg_err_hist = mean(min(err_hist, [], 2));          % 每个时刻取最小误差
avg_err_weight = mean(min(err_weight, [], 2));
avg_err_cbf = mean(min(err_cbf, [], 2));
avg_err_mvdr = mean(min(err_mvdr, [], 2));
avg_err_music = mean(min(err_music, [], 2));

fprintf('\n平均误差统计:\n');
fprintf('平均声强法: %.2f 度\n', avg_err_Iavg);
fprintf('直方图法: %.2f 度\n', avg_err_hist);
fprintf('加权直方图法: %.2f 度\n', avg_err_weight);    
fprintf('CBF: %.2f 度\n', avg_err_cbf);    
fprintf('MVDR: %.2f 度\n', avg_err_mvdr);    
fprintf('MUSIC: %.2f 度\n', avg_err_music);    

