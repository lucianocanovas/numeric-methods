def power_iteration(A: list[list[float]], x0: list[float], tolerance: float, max_iterations: int = 50):
    x = [float(value) for value in x0]
    if not x:
        return 0.0, []

    scale = max(abs(value) for value in x)
    if scale == 0:
        return 0.0, [0.0] * len(x)

    x = [value / scale for value in x]
    dominant_index = max(range(len(x)), key=lambda index: abs(x[index]))
    if x[dominant_index] < 0:
        x = [-value for value in x]

    eigenvalue = 0.0
    for _ in range(max_iterations):
        x_new_raw = [sum(A[i][j] * x[j] for j in range(len(A))) for i in range(len(A))]
        norm = max(abs(value) for value in x_new_raw)
        if norm == 0:
            return 0.0, [0.0] * len(x)

        x_new = [value / norm for value in x_new_raw]
        dominant_index = max(range(len(x_new)), key=lambda index: abs(x_new[index]))
        if x_new[dominant_index] < 0:
            x_new = [-value for value in x_new]

        denominator = sum(value * value for value in x)
        if denominator == 0:
            return 0.0, [0.0] * len(x)

        eigenvalue = sum(x[i] * x_new_raw[i] for i in range(len(x))) / denominator
        err = max(abs(x_new[i] - x[i]) for i in range(len(x)))
        x = x_new
        if err < tolerance:
            return eigenvalue, x

    return eigenvalue, x