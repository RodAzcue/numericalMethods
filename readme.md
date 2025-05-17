### Documentation for Numerical Methods Code

>>The documentation was written by AI, 
>>the code is complately coded by ROdrigo Azcué

---

#### **Overview**
This code implements various numerical methods for solving mathematical problems encountered in scientific computing and engineering. It provides tools for:
- Root-finding (non-linear equations)
- Solving systems of equations
- Interpolation and regression
- Numerical differentiation and integration
- Solving ordinary differential equations (ODEs)
- Statistical analysis

**Dependencies**:
- `sympy`: Symbolic mathematics
- `numpy`: Numerical operations
- `terminalTrueColors`: Terminal color formatting
- `statistics`: Basic statistical functions

---

### **1. Root-Finding Methods**
#### **Class Hierarchy**:
- **`NumSolver`**: Base class for root-finding algorithms.
  - **Subclasses**:
    - `Bisection`: Implements the bisection method.
    - `FalsePosition`: Implements the false position method.
    - `NewtonRaphson`: Implements the Newton-Raphson method.
    - `FixedPoint`: Implements fixed-point iteration.
    - `Secant`: Implements the secant method.

#### **Usage**:
```python
x = symbols('x')
eq = x**3 - 2*x - 5

# Example: Newton-Raphson
solver = NewtonRaphson(eq, [3.0], maxError=0.0001, solution=True)
```

**Key Methods**:
- `printTable()`: Prints iteration steps and errors.
- `getError()`: Calculates approximate error for stopping criteria.
- `getRow()`: Computes the next iteration value.

---

### **2. Systems of Equations**
#### **Class Hierarchy**:
- **`NumeSolverSE`**: Base class for system solvers.
  - **Subclasses**:
    - `FixedPointSE`: Fixed-point iteration for systems.
    - `NewtonRaphsonSE`: Newton-Raphson for systems.

#### **Usage**:
```python
x1, x2 = symbols('x1 x2')
f1 = x1**2 + x2**2 - 4
f2 = x1*x2 - 1

solver = NewtonRaphsonSE([f1, f2], [1.0, 1.0], maxError=0.0001)
```

**Key Features**:
- Uses symbolic Jacobian matrices for Newton-Raphson.
- Handles both linear and non-linear systems.

---

### **3. Interpolation**
#### **Classes**:
- **`Lagrange`**: Constructs Lagrange interpolating polynomials.
- **`Newton`**: Constructs Newton divided difference polynomials.

#### **Usage**:
```python
x_data = [0, 2, 4]
y_data = [1, 3, 7]
interp = Newton(x_data, y_data)
interp.subs(1.5)  # Evaluate at x=1.5
```

---

### **4. Regression**
#### **Class**:
- **`leastSquares`**: Performs linear/non-linear regression.
  - Supports models: `linear`, `exp`, `pow`, `frac`.

#### **Usage**:
```python
x_data = [1, 2, 3]
y_data = [2.1, 3.9, 5.8]
model = leastSquares(x_data, y_data, model='linear')
print(model.eq)  # Prints the regression equation
```

**Output Metrics**:
- Correlation coefficient (`r`).
- Standard error (`stde`).

---

### **5. Numerical Differentiation**
#### **Classes**:
- **`numericalDerivative`**: For equally spaced data.
- **`unequalNumericalDerivative`**: For unequally spaced data.

#### **Usage**:
```python
x_data = [0, 0.1, 0.2]
y_data = [sin(xi) for xi in x_data]
deriv = numericalDerivative(x_data, y_data, degree=1, quality=1)
```

**Methods**:
- `printTable()`: Dispoints derivative estimates at all points.
- `getRow()`: Computes derivative at a specific point.

---

### **6. Numerical Integration**
#### **Methods**:
- **Rectangle Rule**: `rectangleRule()`
- **Trapezoidal Rule**: `trapezoidalRule()`
- **Simpson's 1/3 Rule**: `simpsons13Rule()`
- **Simpson's 3/8 Rule**: `simpsons38Rule()`
- **Gauss Quadrature**: `gaussQuadrature()`

#### **Usage**:
```python
x = symbols('x')
eq = exp(-x**2)
x_data = [0, 0.5, 1.0, 1.5, 2.0]

result = simpsons13Rule(eq, x_data)  # Simpson's 1/3 Rule
```

**Features**:
- Automatic error calculation against symbolic integrals.
- Input validation for method-specific requirements.

---

### **7. Ordinary Differential Equations (ODEs)**
#### **Class**:
- **`RungeKutta`**: Implements RK methods (Euler, Heun, Midpoint, RK4).

#### **Usage**:
```python
x = symbols('x')
y = symbols('y')
ode = y - x  # dy/dx = y - x
solver = RungeKutta(xData=[0, 0.5, 1.0], y0=1, eq=ode, method="rk4")
```

**Supported Methods**:
- `euler`, `heun`, `midpoint`, `ralston`, `rk4`.

---

### **8. Statistical Analysis**
#### **Function**:
- **`statisticsInfo(data)`**: Computes mean, median, standard deviation, etc.

#### **Usage**:
```python
data = [2.1, 3.5, 4.0, 5.2]
statisticsInfo(data)
```

---

### **Execution**
Run the `main()` function to execute predefined examples:
```python
if __name__ == '__main__':
    main()
```

---

If theres an error, make a pull request :)