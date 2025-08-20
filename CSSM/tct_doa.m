function [doa_est, spectrum1] = tct_doa(x, f_bands, fs, M, d, theta_search, f0_opt)
% TCT算法实现宽带相干信号DOA估计
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
    D = 2;                    % 假设信号源数目为2，实际应用中可通过信息准则估计
    
    % ---------------------- 步骤1：粗估角度（ISSM） ----------------------
    K = length(f_bands);            % 频率点数量
    doa_coarse = [];
    Rx = fft(x').';
    for j = 1:K
        % 对每个频率点进行FFT
        Xj = Rx(:, round(f_bands(j)*size(x,2)/fs) );
        Rxj = Xj*Xj'/size(Xj,2); % 协方差矩阵（式4-23的Rx(fj)）
        sigma2_j = estimate_noise(Rxj, D); % 噪声功率估计（小特征值平均）
        Pj = Rxj - sigma2_j*eye(M); % 去噪协方差（式4-23的P(fj)）
        [doa_j, ~] = MUSIC_DOACore(Pj, M, d, c, f_bands(j), theta_search, D); % 粗估DOA
        doa_coarse = [doa_coarse; doa_j];
    end
    doa_coarse = unique(round(doa_coarse*2)/2); % 去重（精度0.5°）
    K_est = length(doa_coarse); % 实际信号数

    % ---------------------- 步骤2：构造聚焦角度集 ----------------------
    BW = 10; % 3dB波束宽度（可调整）
    theta_grid = [];
    for k = 1:K_est
        theta_k = doa_coarse(k);
        theta_grid = [theta_grid, theta_k-0.25*BW:0.5:theta_k+0.25*BW]; % 聚焦范围
    end
    theta_grid = unique(theta_grid); % 去重
    % 对角度网格进行边界截断，确保所有角度值都落在预设的搜索范围内
    theta_grid = theta_grid(theta_grid>=min(theta_search) & theta_grid<=max(theta_search));
    
    % ---------------------- 步骤3：计算子带信号协方差Rs(fj) ----------------------
    Rs_list = cell(K,1); % 存储各子带Rs(fj)
    P_list = cell(K,1); % 存储各子带Pj
    for j = 1:K
        fj = f_bands(j);
        Xj = Rx(:, round(f_bands(j)*size(x,2)/fs) );
        Rxj = Xj*Xj'/size(Xj,2); % 协方差矩阵（式4-23的Rx(fj)） 
        sigma2_j = estimate_noise(Rxj, D);
        Pj = Rxj - sigma2_j*eye(M); P_list{j} = Pj;
        A_fj = construct_steering_vectors(M, d, c, fj, theta_grid); % 导向矩阵A(fj,θ)
        AHA = A_fj'*A_fj; AHA_inv = pinv(AHA);
        Rs_j = AHA_inv*A_fj'*Pj*A_fj*AHA_inv; % 式4-22：Rs(fj) = [A^H A]^{-1}A^H Pj A [A^H A]^{-1}
        Rs_list{j} = Rs_j;
    end
%    % ---------------------- 步骤4：一维搜索最佳聚焦频率f0----------------------
%     f_candidates = linspace(min(f_bands), max(f_bands), 50); % 候选频率
%     min_error = inf; f0_opt = mean(f_bands);
%     for p = 1:length(f_candidates)
%         f0 = f_candidates(p);
%         A_f0 = construct_steering_vectors(M, d, c, f0, theta_grid); % A(f0,θ)
%         AHA0 = A_f0'*A_f0; % AHA0_inv = inv(AHA0);
%         % 计算Rs(f0) = (1/K)ΣRs(fj)（式4-24）
%         Rs_f0 = zeros(size(Rs_list{1}));
%         for j = 1:K, Rs_f0 = Rs_f0 + Rs_list{j}; end
%         Rs_f0 = Rs_f0/K;
%         % 反推P(f0) = A [A^H A] Rs [A^H A] A^H（从式4-22逆推）
%         P_f0 = A_f0*AHA0*Rs_f0*AHA0*A_f0'; 
%         % 计算奇异值误差（式4-25）
%         [~, sigma0, ~] = svd(P_f0); sigma0 = diag(sigma0(1:D));
%         u = zeros(D,1);
%         for j = 1:K
%             [~, sigmaj, ~] = svd(P_list{j}); sigmaj = diag(sigmaj(1:D));
%             u = u + sigmaj;
%         end
%         u = u/K;
%         error = sum(abs(sigma0 - u)); % 误差和
%         if error < min_error, min_error = error; f0_opt = f0; end
%     end
    
    % ---------------------- 步骤5：计算聚焦矩阵T(fj)并合成协方差（式4-19） ----------------------
    R_focus = zeros(M,M);
    for j = 1:K
        fj = f_bands(j);
        A_fj = construct_steering_vectors(M, d, c, fj, theta_grid);
        A_f0 = construct_steering_vectors(M, d, c, f0_opt, theta_grid);
        AHAj = A_fj'*A_fj; AHAj_inv = pinv(AHAj);
        Tj = A_f0*AHAj_inv*A_fj'; % 聚焦矩阵：T(fj)A(fj)=A(f0)（最小二乘解）
        R_focus = R_focus + Tj*P_list{j}*Tj'; % 聚焦后协方差累加
    end
    R_focus = R_focus/K; % 平均

    % ---------------------- 步骤6：最终DOA估计（MUSIC，图中步骤6） ----------------------
    [doa_est, spectrum1] = MUSIC_DOACore(R_focus, M, d, c, f0_opt, theta_search, D);

end

function sigma2 = estimate_noise(R, K)
    % 估计噪声功率：小特征值平均（图中步骤3）
    [~, D] = eig(R);
    eig_vals = diag(D);
    sigma2 = mean(eig_vals(1:end-K)); % 噪声特征值平均
end

function [Aj] = construct_steering_vectors(M, d, c, fj, theta_range)
% 构造参考频率和当前频率的导向矢量矩阵
    N_theta = length(theta_range);
    Aj = zeros(M, N_theta);  % 当前频率fj的导向矢量矩阵
    for i = 1:N_theta
        Aj(:, i) = steering_vector(M, d, c, fj, theta_range(i));
    end
end

function [doa, spectrum] = MUSIC_DOACore(R, M, d, c, f, theta_search, K)
    % MUSIC核心：协方差矩阵→DOA估计（图中步骤2、6）
    [V, ~] = svd(R);
    Un = V(:, K+1:end); % 噪声子空间矩阵
    % 计算空间谱
    spectrum = zeros(size(theta_search));
    for i = 1:length(theta_search)
        theta = theta_search(i);
        a = steering_vector(M, d, c, f, theta);
        spectrum(i) = abs(1/(a'*(Un*Un')*a)); % MUSIC谱公式
    end
    % 提取谱峰
    [~, locs] = findpeaks(spectrum, 'SortStr', 'descend', 'NPeaks', K);
    doa = theta_search(sort(locs(1:K)));
end
