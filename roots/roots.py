from typing import Callable

def bisection(f: Callable, a: float, b: float, tolerance: float, max_iterations: int=50):
    if f(a) * f(b) >= 0:
        return None
    i = 0
    root = []
    r = (a + b) / 2
    root.append(r)
    while abs(f(r)) > tolerance and i < max_iterations:
        if f(a) * f(r) < 0:
            b = r
        else:
            a = r
        r = (a + b) / 2
        i += 1
        root.append(r)
    return root

def regula_falsi(f: Callable, a: float, b: float, tolerance: float, max_iterations: int=50):
    if f(a) * f(b) >= 0:
        return None
    i = 0
    root = []
    r = (a * f(b) - b * f(a)) / (f(b) - f(a))
    root.append(r)
    while abs(f(r)) > tolerance and i < max_iterations:
        if f(a) * f(r) < 0:
            b = r
        else:
            a = r
        r = (a * f(b) - b * f(a)) / (f(b) - f(a))
        i += 1
        root.append(r)
    return root

def secant(f: Callable, x0: float, x1: float, tolerance: float, max_iterations: int=50):
    i = 0
    root = []
    r = x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0))
    root.append(r)
    while abs(f(r)) > tolerance and i < max_iterations:
        x0 = x1
        x1 = r
        r = x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0))
        i += 1
        root.append(r)
    return root

def newton_raphson(f: Callable, x0: float, tolerance: float, max_iterations: int=50):
    i = 0
    root = []
    r = x0 - f(x0) * 1e-5 / (f(x0 + 1e-5) - f(x0))
    root.append(r)
    while abs(f(r)) > tolerance and i < max_iterations:
        x0 = r
        r = x0 - f(x0) * 1e-5 / (f(x0 + 1e-5) - f(x0))
        i += 1
        root.append(r)
    return root

def fixed_point(f: Callable, g: Callable, x0: float, tolerance: float, max_iterations: int=50):
    i = 0
    root = []
    x = g(x0)
    root.append(x)
    while abs(f(x0)) > tolerance and i < max_iterations:
        x0 = x
        x = g(x0)
        root.append(x)
        i += 1
    return root