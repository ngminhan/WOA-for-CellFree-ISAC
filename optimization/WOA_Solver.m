% WOA_Solver.m
% Hàm thực hiện thuật toán Whale Optimization Algorithm (WOA)

function [X_best, fitness_best] = WOA_Solver(fitness_func, dim, bounds, params, H_comm, sensing_beamsteering, gamma_linear)

    % Lấy tham số WOA từ cấu trúc params
    pop_size = params.woa.pop_size;
    max_iter = params.woa.max_iter;
    b = 1; % Hằng số cho spiral (thường là 1)
    
    % --- 1. KHỞI TẠO QUẦN THỂ ---
    pop = zeros(pop_size, dim); % Ma trận vị trí cá voi
    fitness = zeros(pop_size, 1); % Vector Fitness tương ứng
    
    % Khởi tạo ngẫu nhiên trong giới hạn (bounds)
    for i = 1:pop_size
        for j = 1:dim
            % Khởi tạo ngẫu nhiên giữa lower_bound và upper_bound
            pop(i, j) = bounds(j, 1) + (bounds(j, 2) - bounds(j, 1)) * rand();
        end
        % Đánh giá Fitness ban đầu bằng hàm đã tạo ở Bước 2
        fitness(i) = fitness_func(pop(i, :), params, H_comm, sensing_beamsteering, gamma_linear);
    end
    
    % Tìm cá thể tốt nhất ban đầu (Global Best)
    [fitness_best, best_idx] = min(fitness);
    X_best = pop(best_idx, :);
    
    % --- 2. VÒNG LẶP CHÍNH CỦA THUẬT TOÁN WOA ---
    for t = 1:max_iter
        % Cập nhật tham số 'a' (giảm tuyến tính từ 2 về 0)
        a = 2 - t * (2 / max_iter);
        
        for i = 1:pop_size
            r1 = rand(); r2 = rand(); 
            A = 2 * a * r1 - a;       
            C = 2 * r2;               
            l = -1 + 2 * rand();      % Số ngẫu nhiên [-1, 1]
            p = rand();               % Xác suất chuyển đổi cơ chế
            
            % --- CƠ CHẾ CẬP NHẬT VỊ TRÍ ---
            
            if p < 0.5
                % 1. Tận dụng (Exploitation) HOẶC Khám phá (Exploration)
                if abs(A) < 1 
                    % Bao vây Con mồi (Exploitation) 
                    D = abs(C .* X_best - pop(i, :));
                    new_pos = X_best - A .* D;
                else
                    % Tìm kiếm Con mồi (Exploration)
                    rand_idx = randi(pop_size);
                    X_rand = pop(rand_idx, :);
                    D = abs(C .* X_rand - pop(i, :));
                    new_pos = X_rand - A .* D;
                end
            else
                % 2. Tấn công bằng Búp sóng Bong bóng (Spiral movement - Exploitation)
                D_prime = abs(X_best - pop(i, :)); 
                new_pos = D_prime .* exp(b * l) .* cos(2 * pi * l) + X_best;
            end
            
            % --- ÁP DỤNG RÀNG BUỘC GIỚI HẠN (Bounds Check) ---
            % Đảm bảo vị trí mới vẫn nằm trong giới hạn vật lý
            for j = 1:dim
                new_pos(j) = max(new_pos(j), bounds(j, 1));
                new_pos(j) = min(new_pos(j), bounds(j, 2));
            end
            
            % --- ĐÁNH GIÁ VÀ CẬP NHẬT ---
            pop(i, :) = new_pos;
            current_fitness = fitness_func(pop(i, :), params, H_comm, sensing_beamsteering, gamma_linear);
            fitness(i) = current_fitness;
            
            % Cập nhật Global Best
            if current_fitness < fitness_best
                fitness_best = current_fitness;
                X_best = pop(i, :);
            end
        end
        
        % Hiển thị tiến trình
        disp(['Iteration ' num2str(t) '/' num2str(max_iter) ': Best Fitness = ' num2str(fitness_best)]);
    end
end