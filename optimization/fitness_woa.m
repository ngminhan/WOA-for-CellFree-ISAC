% fitness_woa.m
function fitness = fitness_woa(X_vector, params, H_comm, sensing_beamsteering, gamma_linear)

    U = params.U; T = params.T; M_t = params.M_t; N_t = params.N_t; S = U + T; 
    dim_complex = M_t * N_t * S;
    
    % --- 1. TÁI TẠO BEAMFORMING ---
    real_part = X_vector(1:dim_complex);
    imag_part = X_vector(dim_complex+1 : end);
    F_complex = real_part + 1i * imag_part;
    F_streams = reshape(F_complex, [N_t, M_t, S]); 
    F_streams = permute(F_streams, [3, 2, 1]); 
    F_comm = F_streams(1:U, :, :);
    F_sensing = F_streams(U+1:end, :, :);

    % --- 2. TÍNH MỤC TIÊU GỐC ---
    sigmasq_radar_rcs = params.sigmasq_radar_rcs;
    S_SNR = compute_sensing_SNR(sigmasq_radar_rcs, sensing_beamsteering, F_comm, F_sensing);
    objective_term = -S_SNR; 

    % --- 3. TÍNH VÀ ÁP DỤNG HÀM PHẠT ---
    
    % a. Phạt SINR
    penalty_violation_SINR_sum = 0; 
    sigmasq_ue = params.sigmasq_ue;
    SINR_values = compute_SINR(H_comm, F_comm, F_sensing, sigmasq_ue);
    
    for u = 1:U
        penalty_u = max(0, gamma_linear - SINR_values(u));
        penalty_violation_SINR_sum = penalty_violation_SINR_sum + penalty_u;
    end
    penalty_term_SINR = params.lambda_SINR * penalty_violation_SINR_sum;
    
    % b. Phạt Công suất AP tối đa
    P_max = params.P; 
    power_violation_sum = 0; 
    
    for m_t = 1:M_t
        P_tx_m = sum(abs(F_streams(:, m_t, :)).^2, 'all'); 
        penalty_P_m = max(0, P_tx_m - P_max);
        power_violation_sum = power_violation_sum + penalty_P_m;
    end
    
    penalty_term_Power = params.lambda_Power * power_violation_sum;

    % c. Phạt Sensing SNR tối thiểu (S_SNR >= S_th)
    S_th = params.S_th;
    
    % Tính mức vi phạm: max(0, Ngưỡng - S_SNR thực tế)
    penalty_violation_Sensing_sum = max(0, S_th - S_SNR);

    % Áp dụng hệ số phạt
    penalty_term_Sensing = params.lambda_S * penalty_violation_Sensing_sum;
    
    % --- 4. TÍNH GIÁ TRỊ FITNESS CUỐI CÙNG ---
    fitness = objective_term + penalty_term_SINR + penalty_term_Power + penalty_term_Sensing;
    
    if isnan(fitness) || isinf(fitness)
        fitness = 1e10; 
    end
end