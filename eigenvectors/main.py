from eigenvectors import *

A = [[-10, 4], [-4, 0]]
x0 = [-3, 2]
tolerance = 1e-6
eigenvalue, eigenvector = power_iteration(A, x0, tolerance)
print("Dominant eigenvalue:", eigenvalue)
print("Corresponding eigenvector:", eigenvector)