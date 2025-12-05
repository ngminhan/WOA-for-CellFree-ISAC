% reconstruct_F.m
% Hàm chuyển đổi vector tối ưu X_vector (Real) thành các tensor Beamforming (Complex)

function [F_comm, F_sensing] = reconstruct_F(X_vector, params)

    % Lấy các tham số cơ bản từ cấu trúc params
    U = params.U;
    T = params.T;
    M_t = params.M_t;
    N_t = params.N_t;
    S = U + T; % Tổng số luồng
    
    % --- 1. TÁI TẠO BEAMFORMING (X_vector -> F_streams) ---
    
    % Kích thước tổng số phần tử phức
    dim_complex = M_t * N_t * S;
    
    % Tách nửa đầu là Phần Thực và nửa sau là Phần Ảo
    % Lưu ý: X_vector là đầu ra của WOA, nên nó là hàng vector (1 x Dim)
    real_part = X_vector(1:dim_complex);
    imag_part = X_vector(dim_complex+1 : end);
    
    % Khôi phục các phần tử phức
    F_complex = real_part + 1i * imag_part;
    
    % Định hình lại F_complex thành tensor (Streams x M_t x N_t)
    % Quá trình reshape và permute giống hệt trong fitness_woa.m
    F_streams = reshape(F_complex, [N_t, M_t, S]); 
    F_streams = permute(F_streams, [3, 2, 1]); % S x M_t x N_t 
    
    % Tách thành F_comm và F_sensing
    F_comm = F_streams(1:U, :, :);
    F_sensing = F_streams(U+1:end, :, :);

end