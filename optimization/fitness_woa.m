% fitness_woa.m
% Hàm Mục tiêu CẦN TỐI THIỂU HÓA cho WOA: 
% Fitness = -Sensing_SNR + Penalty_SINR + Penalty_Power

function fitness = fitness_woa(X_vector, params, H_comm, sensing_beamsteering, gamma_linear)

    % Lấy các tham số cơ bản từ cấu trúc params
    U = params.U;
    T = params.T;
    M_t = params.M_t;
    N_t = params.N_t;
    S = U + T; % Tổng số luồng
    
    % --- 1. TÁI TẠO BEAMFORMING (X_vector -> F_streams) ---
    
    dim_complex = M_t * N_t * S;
    
    % Tách nửa đầu là Phần Thực và nửa sau là Phần Ảo
    real_part = X_vector(1:dim_complex);
    imag_part = X_vector(dim_complex+1 : end);
    
    % Khôi phục các phần tử phức (Complex Elements)
    F_complex = real_part + 1i * imag_part;
    
    % Định hình lại F_complex thành ma trận/tensor (Streams x M_t x N_t)
    F_streams = reshape(F_complex, [N_t, M_t, S]); 
    F_streams = permute(F_streams, [3, 2, 1]); % S x M_t x N_t (Streams x APs x Antennas)
    
    % Tách thành F_comm và F_sensing
    F_comm = F_streams(1:U, :, :);
    F_sensing = F_streams(U+1:end, :, :);

    % --- 2. TÍNH MỤC TIÊU GỐC: Sensing SNR ---
    
    sigmasq_radar_rcs = params.sigmasq_radar_rcs;
    % Gọi hàm gốc: compute_sensing_SNR.m
    S_SNR = compute_sensing_SNR(sigmasq_radar_rcs, sensing_beamsteering, F_comm, F_sensing);
    
    % Mục tiêu: Tối đa hóa S_SNR. Chuyển thành Tối thiểu hóa (-S_SNR)
    objective_term = -S_SNR; 

    % --- 3. TÍNH VÀ ÁP DỤNG HÀM PHẠT (Penalty Function) ---
    
    % a. Phạt SINR (SINR_u >= gamma_linear)
    penalty_violation_SINR_sum = 0; % SỬA: Biến tích lũy vi phạm SINR
    sigmasq_ue = params.sigmasq_ue;
    % Gọi hàm gốc: compute_SINR.m
    SINR_values = compute_SINR(H_comm, F_comm, F_sensing, sigmasq_ue);
    
    for u = 1:U
        % Cộng mức vi phạm: max(0, ngưỡng - SINR_thực tế)
        penalty_u = max(0, gamma_linear - SINR_values(u));
        % SỬA: Tích lũy vào biến SINR
        penalty_violation_SINR_sum = penalty_violation_SINR_sum + penalty_u;
    end
    
    penalty_term_SINR = params.lambda_SINR * penalty_violation_SINR_sum;
    
    % b. Phạt Công suất AP tối đa (P_tx_m <= P_max)
    P_max = params.P; 
    
    power_violation_sum = 0; % SỬA: KHỞI TẠO biến tích lũy vi phạm Công suất
    
    for m_t = 1:M_t
        % Tính tổng công suất phát của AP thứ m_t (tổng công suất trên tất cả luồng)
        P_tx_m = sum(abs(F_streams(:, m_t, :)).^2, 'all'); 
        
        % Cộng mức vi phạm: max(0, P_thực tế - P_max)
        penalty_P_m = max(0, P_tx_m - P_max);
        
        % SỬA: Tích lũy vào biến Công suất
        power_violation_sum = power_violation_sum + penalty_P_m;
        
        % LƯU Ý: Không tích lũy vào biến SINR (penalty_violation_SINR_sum) ở đây!
    end
    
    % SỬA: Sử dụng biến power_violation_sum đã được tính toán
    penalty_term_Power = params.lambda_Power * power_violation_sum;
    
    % --- 4. TÍNH GIÁ TRỊ FITNESS CUỐI CÙNG ---
    % Cộng các Term Phạt (đã nhân lambda) vào Mục tiêu gốc
    fitness = objective_term + penalty_term_SINR + penalty_term_Power;
    
    % Xử lý trường hợp lỗi (nếu xảy ra NaN hoặc Inf)
    if isnan(fitness) || isinf(fitness)
        fitness = 1e10; % Gán giá trị rất lớn để WOA loại bỏ giải pháp này
    end
end