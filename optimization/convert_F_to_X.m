% convert_F_to_X.m
% Chuyển đổi tensor Beamforming F_init (complex) sang vector X_vector (real)

function [X_vector] = convert_F_to_X(F_streams, params)

    M_t = params.M_t;
    N_t = params.N_t;
    S = params.U + params.T; 
    
    % Kiểm tra F_streams có đúng định dạng S x M_t x N_t không
    if size(F_streams, 1) ~= S || size(F_streams, 2) ~= M_t || size(F_streams, 3) ~= N_t
        % Nếu F_streams chỉ là F_comm (U x M_t x N_t), cần phải điều chỉnh
        % Giả định: F_streams là F_comm (U x M_t x N_t). Cần thêm luồng cảm biến bằng 0.
        if size(F_streams, 1) == params.U && params.T > 0
            F_sensing_zero = zeros(params.T, M_t, N_t);
            F_streams = cat(1, F_streams, F_sensing_zero);
        else
            error('Lỗi: Kích thước F_streams không khớp với hệ thống (S x M_t x N_t).');
        end
    end

    % Làm phẳng F_streams thành vector 1D phức
    F_complex = reshape(permute(F_streams, [3, 2, 1]), [1, M_t * N_t * S]);
    
    dim_complex = M_t * N_t * S;
    
    % Tách phần thực và ảo và ghép lại thành X_vector
    real_part = real(F_complex);
    imag_part = imag(F_complex);
    
    X_vector = [real_part, imag_part]; 
end