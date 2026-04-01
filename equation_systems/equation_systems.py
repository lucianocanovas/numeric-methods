def jacobi(A: list[list[float]], b: list[float], tolerance: float, max_iterations: int = 50):
    x = [0.0] * len(b)
    err = float('inf')
    while err > tolerance:
        x_new = [0.0] * len(x)
        for i in range(len(A)):
            sum_ax = sum(A[i][j] * x[j] for j in range(len(A)) if j != i)
            x_new[i] = (b[i] - sum_ax) / A[i][i]
        err = max(abs(x_new[i] - x[i]) for i in range(len(x)))
        x = x_new
    return x

def alt_jacobi(A: list[list[float]], b: list[float], tolerance: float, max_iterations: int = 50):
    x = [0.0] * len(b)
    D = [[A[i][i] if i == j else 0.0 for j in range(len(A))] for i in range(len(A))]
    B = [[A[i][j] - D[i][j] for j in range(len(A))] for i in range(len(A))]
    T = [[-B[i][j] / D[i][i] if D[i][i] != 0 else 0.0 for j in range(len(A))] for i in range(len(A))]
    err = float('inf')
    while err > tolerance:
        x_new = [0.0] * len(x)
        for i in range(len(A)):
            sum_tx = sum(T[i][j] * x[j] for j in range(len(A)))
            x_new[i] = sum_tx + b[i] / D[i][i]
        err = max(abs(x_new[i] - x[i]) for i in range(len(x)))
        x = x_new
    return x

def gauss_seidel(A: list[list[float]], b: list[float], tolerance: float, max_iterations: int = 50):
    x = [0.0] * len(b)
    err = float('inf')
    while err > tolerance:
        x_new = x.copy()
        for i in range(len(A)):
            sum_ax = sum(A[i][j] * x_new[j] for j in range(len(A)) if j != i)
            x_new[i] = (b[i] - sum_ax) / A[i][i]
        err = max(abs(x_new[i] - x[i]) for i in range(len(x)))
        x = x_new
    return x