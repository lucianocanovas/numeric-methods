function eigenvectors()
    format long;
    A = [-10, 4; -4, 0];
    x0 = [-3; 2];
    tolerance = 1e-4;
    max_iterations = 50;
    lambda = power_iteration(A, x0, tolerance, max_iterations);
end

function [lambda, x] = power_iteration(A, x0, tolerance, max_iterations)
    x = x0 / norm(x0, inf);
    [~, index] = max(abs(x));
    if x(index) < 0
        x = -x;
    endif

    lambda = 0;
    for i = 1:max_iterations
        A_new = A * x;
        norm_new = norm(A_new, inf);
        if norm_new == 0
            lambda = 0;
            x = zeros(size(x));
            return;
        endif

        x_new = A_new / norm_new;
        [~, index] = max(abs(x_new));
        if x_new(index) < 0
            x_new = -x_new;
        endif

        lambda_new = (x' * A_new) / (x' * x);
        if norm(x_new - x, inf) < tolerance
            lambda = lambda_new
            x = x_new
            return;
        endif

        x = x_new;
        lambda = lambda_new;
    endfor
end
