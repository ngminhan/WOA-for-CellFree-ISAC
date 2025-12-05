cvx_path = '/Users/nguyenminhan/Desktop/Documents/GitHub/cvx'; % Path of the CVX installation
addpath(genpath(cvx_path)); % Add CVX to the path
addpath(genpath('./')); % Add paths of subfolders

%% Simulation for Fig. 3
clear all;
eval('sim_params'); % load the parameters

% --- THÊM: Tính toán Dim và Bounds cho WOA (Cần thiết sau clear all) ---
S = params.U + params.T; % Tổng số luồng Beamforming 
Dim = 2 * S * params.M_t * params.N_t; 
max_amp = sqrt(params.M_t * params.P); 
lower_bound = -max_amp; upper_bound = max_amp; 
bounds = repmat([lower_bound, upper_bound], Dim, 1);
% ---------------------------------------------------------------------

output_file = 'power_data';
% THAY ĐỔI: Thêm Dim và bounds vào lệnh gọi
simulation(params, output_file, Dim, bounds); 
plot_results_power;


%% Simulation for Fig. 4
clear all;
eval('sim_params'); % load the parameters
params.P_comm_ratio = [0.5]; % Fixed power ratio for communications

% --- THÊM: Tính toán Dim và Bounds cho WOA (Đã có logic, chỉ cần giữ) ---
S = params.U + params.T; % Tổng số luồng Beamforming 
Dim = 2 * S * params.M_t * params.N_t; 
% Định nghĩa Giới hạn Tìm kiếm (Bounds) 
max_amp = sqrt(params.M_t * params.P); 
lower_bound = -max_amp; upper_bound = max_amp; 
% Tạo ma trận Bounds [Dim x 2] 
bounds = repmat([lower_bound, upper_bound], Dim, 1);
% ---------------------------------------------------------------------

% Generate samples with a pre-determined range of
% minimum distance between the target and closest UE
for min_dist = 0:5:45   
    params.geo.min_dist = min_dist;
    params.geo.max_dist = min_dist+5;
    
    output_file = strcat('dist_data', num2str(min_dist));
    % THAY ĐỔI: Thêm Dim và bounds vào lệnh gọi
    simulation(params, output_file, Dim, bounds); 
end
plot_results_dist;