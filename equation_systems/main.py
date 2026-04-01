from equation_systems import *

A = [[50, -25, 0, 0], [-25, 50, -25, 0], [0, -25, 50, -25], [0, 0, -25, 50]]
b = [10, 20, 20, 10]
tolerance = 1e-6
jacobi_solution = jacobi(A, b, tolerance)
alt_jacobi_solution = alt_jacobi(A, b, tolerance)
gauss_seidel_solution = gauss_seidel(A, b, tolerance)
print("Jacobi solution:", jacobi_solution)
print("Alternative Jacobi solution:", alt_jacobi_solution)
print("Gauss-Seidel solution:", gauss_seidel_solution)