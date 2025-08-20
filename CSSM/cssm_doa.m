function [doa_est, spectrum1, spectrum2] = cssm_doa(x, f, fs, M, d, theta_range, f0_opt)
% CSSM算法实现宽带相干信号DOA估计
% 输入参数:
%   x: 接收信号矩阵，大小为M×N，M为阵元数，N为采样点数
%   f: 频率点向量
%   fs: 采样频率
%   M: 阵元数量
%   d: 阵元间距
%   theta_range: 角度搜索范围
% 输出参数:
%   doa_est: 估计的波达方向
%   spectrum: 空间谱
    
    % 参数设置
    c = 1500;                 % 声速(m/s)，若是电磁波可改为光速
    f0 = f0_opt;             % 参考频率选为中心频率
    D = 2;                    % 假设信号源数目为2，实际应用中可通过信息准则估计
    
    % 1. 计算每个频率点的协方差矩阵
    K = length(f);            % 频率点数量
    R = cell(K, 1);           % 存储各频率点的协方差矩阵
    
    for j = 1:K
        % 对每个频率点进行FFT
        Xj = fft(x, [], 2);
        Xj = Xj(:, round(f(j)*size(x,2)/fs) );
        
        % 估计协方差矩阵
        R{j} = Xj * Xj' / size(Xj, 2);
    end
    
    % 2. 构造聚焦矩阵
    T = cell(K, 1);           % 存储各频率点的聚焦矩阵
    
    for j = 1:K
        % 构造参考频率和当前频率的导向矢量矩阵
        [A0, Aj] = construct_steering_vectors(M, d, c, f0, f(j), theta_range);
        
        % 基于最小二乘法构造聚焦矩阵 T = A0 * Aj^+
        T{j} = A0 * pinv(Aj);
    end
    
    % 3. 计算聚焦后的协方差矩阵
    Ry = zeros(M, M);
    for j = 1:K
        Ry = Ry + T{j} * R{j} * T{j}';
    end
    Ry = Ry / K;
    
    % 4. 特征分解获取信号和噪声子空间
%     [V, D_eig] = eig(Ry);
%     [~, idx] = sort(diag(D_eig), 'descend');
%     V = V(:, idx);
    [V, ~] = svd(Ry);
    
    % 分离信号子空间和噪声子空间
    Us = V(:, 1:D);           % 信号子空间
    Un = V(:, D+1:end);       % 噪声子空间
    
    % 5. 计算MUSIC空间谱
    spectrum1 = zeros(size(theta_range));
    spectrum2 = zeros(size(theta_range));
    for i = 1:length(theta_range)
        theta = theta_range(i);
        a = steering_vector(M, d, c, f0, theta);  % 参考频率下的导向矢量
        spectrum1(i) = abs(1 / (a' * Un * Un' * a));
        spectrum2(i) = abs((a' * Us * Us' * a));
    end
    
    % 6. 寻找谱峰，估计DOA
    [~, peaks_idx] = findpeaks(spectrum1, 'SortStr', 'descend', 'NPeaks', D);
    doa_est = theta_range(sort(peaks_idx(1:D)));
end

function [A0, Aj] = construct_steering_vectors(M, d, c, f0, fj, theta_range)
% 构造参考频率和当前频率的导向矢量矩阵
    N_theta = length(theta_range);
    A0 = zeros(M, N_theta);  % 参考频率f0的导向矢量矩阵
    Aj = zeros(M, N_theta);  % 当前频率fj的导向矢量矩阵
    
    for i = 1:N_theta
        A0(:, i) = steering_vector(M, d, c, f0, theta_range(i));
        Aj(:, i) = steering_vector(M, d, c, fj, theta_range(i));
    end
end

