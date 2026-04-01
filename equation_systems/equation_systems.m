function equation_systems
    A = [50, -25, 0, 0;-25, 50, -25, 0; 0, -25, 50, -25; 0, 0, -25, 50];
    B = [10, 20, 20, 10]';
    alt_x = zeros(4, 1);
    tolerance = 1e-6;
    disp('Jacobi method:');
    Jacobi(A, B, alt_x, tolerance);
    disp('Alternative Jacobi method:');
    alt_jaboci(A, B, alt_x, tolerance);
    disp("Verification:")
    Ax = A * [1.2, 2.0, 2.0, 1.2]';
    disp(Ax);
    A = [-3, 1, -2; 4, -5, 0; 1, -3, 6];
    B = [-2, 5, 6];
    alt_x = zeros(3, 1);
    Jacobi(A, B, alt_x, tolerance)
end

function Jacobi(A, b, x, tolerance)
    err = Inf;
    while err > tolerance
        x_new = zeros(size(x));
        for i = 1:length(x)
            sum = A(i, :) * x - A(i, i) * x(i);
            x_new(i) = (b(i) - sum) / A(i, i);
        end
        err = norm(x_new - x, Inf);
        x = x_new;
    end
    disp('Solution:');
    disp(x);
end
    
function alt_jaboci(A, b, x, tolerance)
    D = diag(diag(A));
    B = A - D;
    T = -D\B;
    C = D\b;
    err = Inf;
    while err > tolerance
        x_new = T * x + C;
        err = norm(x_new - x, Inf);
        x = x_new;
    end
    disp('Solution:');
    disp(x);
end

function gauss_seidel(A, b, x, tolerance)
    err = Inf;
    while err > tolerance
        x_old = x;
        for i = 1:length(x)
            sum = A(i, :) * x - A(i, i) * x(i);
            x(i) = (b(i) - sum) / A(i, i);
        end
        err = norm(x - x_old, Inf);
    end
    disp('Solution:');
    disp(x);
end