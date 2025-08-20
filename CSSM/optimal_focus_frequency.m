function [f0_opt, S_values] = optimal_focus_frequency(M, d, c, f_candidates, f_range, K, theta_grid)
% 计算最佳聚焦频率
    % 输入参数:
    %   M: 阵元数量
    %   d: 阵元间距
    %   c: 波速
    %   f_candidates: 候选聚焦频率集合
    %   f_range: 信号可能的频率范围
    %   K: 信号源数目
    %   theta_grid: 角度网格，用于构造方向矩阵
    % 输出参数:
    %   f0_opt: 最佳聚焦频率
    %   S_values: 各候选频率对应的加权和
    
    % 步骤1: 预处理计算u_i
    J = length(f_range);  % 频率点数量
    u = zeros(K, 1);      % 存储累积奇异值
    
    % 对每个频率点计算奇异值并累积
    for j = 1:J
        fj = f_range(j);
        % 构造该频率下的方向矩阵A(fj)
        A_fj = construct_direction_matrix(M, d, c, fj, theta_grid);
        % 计算奇异值
        [~, sigma, ~] = svd(A_fj');
        sigma = diag(sigma);
        % 取前K个较大的奇异值并累加到u中
        for i = 1:K
            if i <= length(sigma)
                u(i) = u(i) + sigma(i);
            end
        end
    end
    
    % 步骤2: 对每个候选聚焦频率计算加权和S(f0)
    num_candidates = length(f_candidates);
    S_values = zeros(num_candidates, 1);
    
    for p = 1:num_candidates
        f0 = f_candidates(p);
        % 构造参考频率下的方向矩阵A(f0)
        A_f0 = construct_direction_matrix(M, d, c, f0, theta_grid);
        % 计算奇异值
        [~, sigma0, ~] = svd(A_f0);
        sigma0 = diag(sigma0);
        
        % 计算加权和S(f0)
        S = 0;
        for i = 1:K
            if i <= length(sigma0)
                S = S + sigma0(i) * u(i);
            end
        end
        S_values(p) = S;
    end
    
    % 步骤3: 找到使S(f0)最大的最佳聚焦频率
    [~, idx_opt] = max(S_values);
    f0_opt = f_candidates(idx_opt);
end

function A = construct_direction_matrix(M, d, c, f, theta_grid)
% 构造方向矩阵A(f)
    % M: 阵元数量
    % d: 阵元间距
    % c: 波速
    % f: 频率
    % theta_grid: 角度网格
    N_theta = length(theta_grid);
    A = zeros(M, N_theta);
    
    for i = 1:N_theta
        theta = theta_grid(i);
        % 构造导向矢量
        theta_rad = theta * pi / 180;
        A(:, i) = exp(1i * 2 * pi * f * d * (0:M-1)' * sin(theta_rad) / c);
    end
end
