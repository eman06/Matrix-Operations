# Matrix Operations

This C++ project features a `Matrix` class designed for handling 2D matrices of floating-point numbers. The class supports various operations and custom memory management, including matrix addition, subtraction, multiplication, and more.

## Features

- **Custom Memory Management**:
  - Implement custom memory allocation and deallocation.
  - Overload `new` and `delete` operators.

- **Matrix Operations**:
  - **Arithmetic Operators**: Overload `+`, `-`, `*` for matrix operations.
  - **Scalar Multiplication**: Implemented using the `*` operator.
  - **Matrix Inversion**: Using Gauss-Jordan elimination.
  - **Determinant Calculation**: Using LU decomposition.
  - **Input/Output Operators**: Overload `>>`, `<<` for matrix I/O.
  - **Comparison Operators**: Overload `==`, `!=` for equality checks.
  - **Element Access**: Overload `()` for accessing and modifying elements with boundary checking.
  - **Unary Negation**: Overload `-` to negate all matrix elements.

- **Additional Features**:
  - Implement functions for matrix transpose, inverse, and eigenvalues.
  
- **File Handling**:
  - Store matrix results in a text file.
  - Reload matrices from the text file upon program restart.

## Sample Usage

```cpp
Matrix m1(3, 3); // Create a 3x3 matrix
Matrix m2(3, 3); // Create another 3x3 matrix

cin >> m1; // Input matrix
m2 = m1 + m2; // Add matrices
cout << m2; // Output matrix

Matrix m3 = m1 * 2.0; // Scalar multiplication
Matrix m4 = m1 * m2; // Matrix multiplication

Matrix m5 = m1.inverse(); // Matrix inversion
double det = m1.determinant(); // Determinant calculation

// Save matrix to file
m1.saveToFile("matrix.txt");

// Load matrix from file
Matrix m6;
m6.loadFromFile("matrix.txt");
