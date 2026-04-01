function roots
    format long;
    f = @(x) x^2 -3;
    g = @(x) (-1/3)*x^2 + x + 1;
    disp("Bisection method:");
    bisection_root = bisection(f, 0, 2, 0.0001);
    disp("Regula Falsi method:");
    regula_falsi_root = regula_falsi(f, 0, 2, 0.0001);
    disp("Secant method:");
    secant_root = secant(f, 0, 2, 0.0001);
    disp("Newton-Raphson method:");
    newton_raphson_root = newton_raphson(f, 0, 0.0001);
    disp("Fixed-point iteration method:");
    fixed_point_root = fixed_point(f, g, 2, 0.0001); 
end

function root = bisection(f, a, b, tolerance)
    if f(a) * f(b) >= 0
        disp('The function must have different signs at the endpoints.');
        root = [];
    end
    r = (a + b)/2
    root = [r];
    while abs(f(r)) > tolerance
        if f(a) * f(r) < 0
            b = r;
        endif
        if f(b) * f(r) < 0
            a = r;
        endif
        r = (a + b)/2
        root = [root, r];
    endwhile
end

function root = regula_falsi(f, a, b, tolerance)
    if f(a) * f(b) >= 0
        disp('The function must have different signs at the endpoints.');
        root = [];
    end
    r = (a * f(b) - b * f(a)) / (f(b) - f(a))
    root = [r];
    while abs(f(r)) > tolerance
        if f(a) * f(r) < 0
            b = r;
        endif
        if f(b) * f(r) < 0
            a = r;
        endif
        r = (a * f(b) - b * f(a)) / (f(b) - f(a))
        root = [root, r];
    endwhile
end

function root = secant(f, x0, x1, tolerance)
    r = x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0))
    root = [r];
    while abs(f(r)) > tolerance
        x0 = x1;
        x1 = r;
        r = x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0))
        root = [root, r];
    endwhile
end

function root = newton_raphson(f, x0, tolerance)
    r = x0 - f(x0) * 1e-5 / (f(x0 + 1e-5) - f(x0))
    root = [r];
    while abs(f(r)) > tolerance
        x0 = r;
        r = x0 - f(x0) * 1e-5 / (f(x0 + 1e-5) - f(x0))
        root = [root, r];
    endwhile
end

function root = fixed_point(f, g, x0, tolerance)
    x = g(x0)
    root = [x];
    while abs(f(x0)) > tolerance
        x0 = x;
        x = g(x0)
        root = [root, x];
    endwhile
end