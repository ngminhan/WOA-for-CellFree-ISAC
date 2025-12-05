%% Simulation
function results = simulation(params, output_filename, Dim, bounds)

    save_filename = output_filename;
    results = cell(1, params.repetitions);

    % Chuyển đổi ngưỡng SINR (gamma) từ dB sang tỉ lệ tuyến tính 
    target_SINR_dB = 5; 
    % Giả sử ngưỡng mục tiêu là 5 dB, bạn có thể điều chỉnh 
    target_gamma_linear = 10^(target_SINR_dB / 10);
    
    % Parfor Supported loop
    parfor rep=1:params.repetitions
        fprintf('\n Repetition %i:', rep)
        [UE_pos, AP_pos, target_pos] = generate_positions(params.T, ...
            params.U, params.M_t, params.geo.line_length, ...
            params.geo.target_y, params.geo.UE_y, ...
            params.geo.min_dist, params.geo.max_dist);
        
        results{rep}.P_comm_ratio = params.P_comm_ratio;
        results{rep}.AP = AP_pos;
        results{rep}.UE = UE_pos;
        results{rep}.Target = target_pos;
    
        % Channel generation
        H_comm = LOS_channel(AP_pos, UE_pos, params.N_t);
        % H_sensing = radar_LOS_channel(target_pos, AP_pos, AP_pos, N_t, N_t, params.sigmasq_radar_rcs);
            
        % results{rep}.H_comm = H_comm;
        % results{rep}.H_sensing = H_sensing;
        
        %% Compute sensing beams
        [sensing_angle, ~] = compute_angle_dist(AP_pos, target_pos);
        sensing_beamsteering = beamsteering(sensing_angle.', params.N_t); % (Target x M_t x N_t)
        
        % Conjugate BF - Power normalized beamsteering
        F_sensing_CB_norm = sensing_beamsteering*sqrt(1/params.N_t);
        
        % Nullspace BF
        F_sensing_NS_norm = beam_nulling(H_comm, sensing_beamsteering);
    
        for p_i = 1:length(params.P_comm_ratio)
    
            % Power ratio of comm and sensing
            P_comm = params.P * params.P_comm_ratio(p_i);
            P_sensing = params.P * (1-params.P_comm_ratio(p_i));
    
            F_sensing_CB = F_sensing_CB_norm * sqrt(P_sensing);
            F_sensing_NS = F_sensing_NS_norm * sqrt(P_sensing);
            
            solution_counter = 1;
    
    
            
            %% NS Sensing - RZF Comm (Baseline)
            F_star_RZF = beam_regularized_zeroforcing(H_comm, P_comm, params.sigmasq_ue)*sqrt(P_comm);
            results{rep}.power{p_i}{solution_counter} = compute_metrics(H_comm, F_star_RZF, params.sigmasq_ue, sensing_beamsteering, F_sensing_NS, params.sigmasq_radar_rcs);
            results{rep}.power{p_i}{solution_counter}.name = 'NS+RZF';
            solution_counter = solution_counter + 1;
    
            %% NS Sensing - Opt Comm (Baseline)
            wrapped_objective = @(gamma) opt_comm_SOCP_vec(H_comm, params.sigmasq_ue, P_comm, F_sensing_NS, gamma);
            [F_star_SOCP_NS, SINR_min_SOCP_NS] = bisection_SINR(params.bisect.low, params.bisect.high, params.bisect.tol, wrapped_objective);
            results{rep}.power{p_i}{solution_counter} = compute_metrics(H_comm, F_star_SOCP_NS, params.sigmasq_ue, sensing_beamsteering, F_sensing_NS, params.sigmasq_radar_rcs);
            results{rep}.power{p_i}{solution_counter}.name = 'NS+OPT';
            results{rep}.power{p_i}{solution_counter}.min_SINR_opt = SINR_min_SOCP_NS;
            solution_counter = solution_counter + 1;
            
            %% CB Sensing - OPT Comm (Baseline)
            wrapped_objective = @(gamma) opt_comm_SOCP_vec(H_comm, params.sigmasq_ue, P_comm, F_sensing_CB, gamma);
            [F_star_SOCP_CB, SINR_min_SOCP_CB] = bisection_SINR(params.bisect.low, params.bisect.high, params.bisect.tol, wrapped_objective);
            results{rep}.power{p_i}{solution_counter} = compute_metrics(H_comm, F_star_SOCP_CB, params.sigmasq_ue, sensing_beamsteering, F_sensing_CB, params.sigmasq_radar_rcs);
            results{rep}.power{p_i}{solution_counter}.name = 'CB+OPT';
            results{rep}.power{p_i}{solution_counter}.min_SINR_opt = SINR_min_SOCP_CB;
            solution_counter = solution_counter + 1;
            
            %% JSC (WOA Solver)
            % Sử dụng SINR_min_SOCP_NS làm ngưỡng SINR tối thiểu (gamma) cho WOA
            sens_streams = 1;
            
            % GỌI WOA SOLVER
            % Hàm @fitness_woa cần được truyền tất cả các tham số cấu trúc
            [X_best_woa, ~] = WOA_Solver(@fitness_woa, Dim, bounds, params, H_comm, sensing_beamsteering, SINR_min_SOCP_NS);
            
            % TÁI TẠO BEAMFORMING VECTORS
            [F_comm_woa, F_sensing_woa] = reconstruct_F(X_best_woa, params);

            % LƯU KẾT QUẢ WOA VÀO CẤU TRÚC KẾT QUẢ
            results{rep}.power{p_i}{solution_counter} = compute_metrics(H_comm, F_comm_woa, params.sigmasq_ue, sensing_beamsteering, F_sensing_woa, params.sigmasq_radar_rcs);
            results{rep}.power{p_i}{solution_counter}.feasible = true; % WOA luôn trả về kết quả
            results{rep}.power{p_i}{solution_counter}.name = strcat('JSC+WOA_',num2str(sens_streams));
            % SSNR_opt là giá trị Sensing SNR đạt được sau cùng
            results{rep}.power{p_i}{solution_counter}.SSNR_opt = results{rep}.power{p_i}{solution_counter}.SSNR; 
            solution_counter = solution_counter + 1;
        end
    end
    
    %% Save Results
    output_folder = './output/';
    if ~exist(output_folder)
        mkdir(output_folder);
    end
    save(strcat(output_folder, save_filename, '.mat'));
end