import matplotlib.pyplot as plt

from roots import *

f = lambda x: x - 2**(-x)
tolerance = 0.0001
bisection_root = bisection(f, 0, 1, tolerance)
regula_falsi_root = regula_falsi(f, 0, 1, tolerance)
secant_root = secant(f, 0, 1, tolerance)
newton_raphson_root = newton_raphson(f, 0, tolerance)
fixed_point_root = fixed_point(f, lambda x: 2**(-x), 0, tolerance)

plt.figure(figsize=(8, 5))
plt.plot(range(len(bisection_root)), bisection_root, label="Bisection")
plt.plot(range(len(regula_falsi_root)), regula_falsi_root, label="Regula Falsi")
plt.plot(range(len(secant_root)), secant_root, label="Secant")
plt.plot(range(len(newton_raphson_root)), newton_raphson_root, label="Newton-Raphson")
plt.plot(range(len(fixed_point_root)), fixed_point_root, label="Fixed-point")
plt.xlabel("Number of iterations")
plt.ylabel("Root approximation")
plt.legend()
plt.tight_layout()
plt.savefig("roots/roots_convergence.png", dpi=150)