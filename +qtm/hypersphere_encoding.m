
function hypersphere_encoding()

% https://chatgpt.com/share/67c68f35-fe90-8005-8c81-c981fb5e2c3a

    n = 3; % Dimension of the hypersphere
    point = randn(n, 1); % Generate a random point
    point = point / norm(point); % Normalize to lie on the unit sphere

    binary_cartesian = cartesian_to_binary(point, 8);
    binary_spherical = spherical_to_binary(point, 8);

    disp('Point:');
    disp(point');
    disp('Binary Cartesian Encoding:');
    disp(binary_cartesian);
    disp('Binary Spherical Encoding:');
    disp(binary_spherical);
end

function bin_str = float_to_binary(value, num_bits, lower_bound, upper_bound)
    % Quantize a float to a fixed number of bits
    scaled_value = round(((value - lower_bound) / (upper_bound - lower_bound)) * (2^num_bits - 1));
    bin_str = dec2bin(scaled_value, num_bits);
end

function bin_str = cartesian_to_binary(point, num_bits)
    % Convert Cartesian coordinates of a hypersphere point to a Boolean string
    bin_str = '';
    for i = 1:length(point)
        bin_str = strcat(bin_str, float_to_binary(point(i), num_bits, -1, 1));
    end
end

function angles = cartesian_to_spherical(point)
    % Convert Cartesian coordinates to spherical coordinates
    n = length(point);
    angles = zeros(1, n - 1);
    
    for i = 1:(n - 1)
        norm_val = norm(point(i:end));
        if norm_val ~= 0
            angles(i) = acos(point(i) / norm_val);
        else
            angles(i) = 0;
        end
    end
end

function bin_str = spherical_to_binary(point, num_bits)
    % Convert a point's spherical coordinates to a Boolean string
    angles = cartesian_to_spherical(point);
    bin_str = '';
    for i = 1:length(angles)
        bin_str = strcat(bin_str, float_to_binary(angles(i), num_bits, 0, pi));
    end
end
