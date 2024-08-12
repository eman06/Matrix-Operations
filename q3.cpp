#include <iostream>
#include <string>
#include <cmath>
#include <cstring>  // For memcpy

using namespace std;

class Matrix {
private:
    int row;
    int col;
    float** mat;

public:
    // Constructors
    Matrix() : row(0), col(0), mat(nullptr) {}

    Matrix(int r, int c) : row(r), col(c) {
        mat = new float* [row];
        for (int i = 0; i < row; ++i) {
            mat[i] = new float[col]();
        }
    }

    Matrix(const Matrix& other) : row(other.row), col(other.col) {
        mat = new float* [row];
        for (int i = 0; i < row; ++i) {
            mat[i] = new float[col];
            for (int j = 0; j < col; ++j) {
                mat[i][j] = other.mat[i][j];
            }
        }
    }

    ~Matrix() {
        for (int i = 0; i < row; ++i) {
            delete[] mat[i];
        }
        delete[] mat;
    }

    // Getters
    int getrow() const { return row; }
    int getcol() const { return col; }
    float** getmat() const { return mat; }

    // Matrix Operations
    Matrix operator+(const Matrix& other) const {
        if (this->row != other.row || this->col != other.col) {
            cout << "Matrix dimensions must be the same for addition." << endl;
            return Matrix(0, 0); // Return an empty matrix to indicate failure
        }

        Matrix sum(this->row, this->col);
        for (int i = 0; i < this->row; ++i) {
            for (int j = 0; j < this->col; ++j) {
                sum.mat[i][j] = this->mat[i][j] + other.mat[i][j];
            }
        }
        return sum;
    }

    Matrix operator-(const Matrix& other) const {
        if (this->row != other.row || this->col != other.col) {
            cout << "Matrix dimensions must be the same for subtraction." << endl;
            return Matrix(0, 0); // Return an empty matrix to indicate failure
        }

        Matrix diff(this->row, this->col);
        for (int i = 0; i < this->row; ++i) {
            for (int j = 0; j < this->col; ++j) {
                diff.mat[i][j] = this->mat[i][j] - other.mat[i][j];
            }
        }
        return diff;
    }

    Matrix operator*(float scalar) const {
        Matrix m(row, col);
        for (int i = 0; i < this->row; ++i) {
            for (int j = 0; j < this->col; ++j) {
                m.mat[i][j] = this->mat[i][j] * scalar;
            }
        }
        return m;
    }

    Matrix operator*(const Matrix& other) const {
        if (this->col != other.row) {
            cout << "Matrices cannot be multiplied due to incompatible dimensions." << endl;
            return Matrix(0, 0); // Return an empty matrix to indicate failure
        }

        Matrix result(this->row, other.col);
        for (int i = 0; i < this->row; ++i) {
            for (int j = 0; j < other.col; ++j) {
                result.mat[i][j] = 0;
                for (int k = 0; k < this->col; ++k) {
                    result.mat[i][j] += this->mat[i][k] * other.mat[k][j];
                }
            }
        }
        return result;
    }

    Matrix transpose() const {
        Matrix m(col, row);
        for (int i = 0; i < this->row; ++i) {
            for (int j = 0; j < this->col; ++j) {
                m.mat[j][i] = this->mat[i][j];
            }
        }
        return m;
    }

    Matrix inverse() const {
        if (row != col) {
            cout << "\nError: Matrix must be square to find its inverse!" << endl;
            return Matrix(); // Return an empty matrix to indicate failure
        }

        int n = row;
        Matrix result(n, n);

        // Create augmented matrix [A|I]
        float** I = new float* [n];
        for (int i = 0; i < n; ++i) {
            I[i] = new float[2 * n]();
            for (int j = 0; j < n; ++j) {
                I[i][j] = mat[i][j];
                I[i][j + n] = (i == j) ? 1 : 0;
            }
        }

        // Apply Gaussian elimination
        for (int i = 0; i < n; ++i) {
            // Find pivot
            int maxRow = i;
            for (int k = i + 1; k < n; ++k) {
                if (fabs(I[k][i]) > fabs(I[maxRow][i])) {
                    maxRow = k;
                }
            }

            if (fabs(I[maxRow][i]) < 1e-9) {
                cout << "\nMatrix is singular!" << endl;
                return Matrix(); // Return an empty matrix to indicate failure
            }

            // Swap rows
            if (i != maxRow) {
                swap(I[i], I[maxRow]);
            }

            // Normalize pivot row
            float pivot = I[i][i];
            for (int j = 0; j < 2 * n; ++j) {
                I[i][j] /= pivot;
            }

            // Eliminate other rows
            for (int k = 0; k < n; ++k) {
                if (k != i) {
                    float factor = I[k][i];
                    for (int j = 0; j < 2 * n; ++j) {
                        I[k][j] -= factor * I[i][j];
                    }
                }
            }
        }

        // Extract the inverse matrix
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                result.mat[i][j] = I[i][j + n];
            }
        }

        // Clean up
        for (int i = 0; i < n; ++i) {
            delete[] I[i];
        }
        delete[] I;

        return result;
    }

    float determinant() const {
        if (row != col) {
            cout << "Matrix must be square to compute the determinant." << endl;
            return 0;
        }
        float** matrix = new float* [row];
        for (int i = 0; i < row; ++i) {
            matrix[i] = new float[col];
            for (int j = 0; j < col; ++j) {
                matrix[i][j] = mat[i][j];
            }
        }
        float det = determinant(matrix, row);

        // Clean up
        for (int i = 0; i < row; ++i) {
            delete[] matrix[i];
        }
        delete[] matrix;

        return det;
    }
    void decomposeQR(float A[10][10], float Q[10][10], float R[10][10]) {
        int n = row;

        // Initialize Q and R matrices
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                Q[i][j] = 0;
                R[i][j] = 0;
            }
        }

        // QR decomposition using Gram-Schmidt process
        for (int i = 0; i < n; ++i) {
            float norm = 0;
            for (int k = 0; k < n; ++k) {
                norm += A[k][i] * A[k][i];
            }
            norm = sqrt(norm);

            for (int k = 0; k < n; ++k) {
                Q[k][i] = A[k][i] / norm;
            }
            R[i][i] = norm;

            for (int j = i + 1; j < n; ++j) {
                float r = 0;
                for (int k = 0; k < n; ++k) {
                    r += Q[k][i] * A[k][j];
                }
                R[i][j] = r;
                for (int k = 0; k < n; ++k) {
                    A[k][j] -= r * Q[k][i];
                }
            }
        }
    }

    float* eigenvalues(int maxIterations = 1000, float tolerance = 1e-6) {
        int n = row;
        float A[10][10];
        memcpy(A, mat, sizeof(A));
        float Q[10][10], R[10][10];

        // Iterate to converge to the eigenvalues
        for (int iter = 0; iter < maxIterations; ++iter) {
            decomposeQR(A, Q, R);

            float newA[10][10] = { 0 };
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    for (int k = 0; k < n; ++k) {
                        newA[i][j] += R[i][k] * Q[k][j];
                    }
                }
            }

            // Check convergence
            bool converged = true;
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    if (fabs(newA[i][j] - A[i][j]) > tolerance) {
                        converged = false;
                        break;
                    }
                }
                if (!converged) break;
            }

            memcpy(A, newA, sizeof(newA));
            if (converged) break;
        }

        float* eigenvals = new float[n];
        for (int i = 0; i < n; ++i) {
            eigenvals[i] = A[i][i];
        }

        return eigenvals;
    }

private:
    void getMinor(float** matrix, float** min, int rrow, int rcol, int n) const {
        int minorRow = 0;
        for (int i = 0; i < n; ++i) {
            if (i == rrow) continue;
            int minorCol = 0;
            for (int j = 0; j < n; ++j) {
                if (j == rcol) continue;
                min[minorRow][minorCol++] = matrix[i][j];
            }
            ++minorRow;
        }
    }

    float determinant(float** matrix, int n) const {
        if (n == 1) {
            return matrix[0][0];
        }
        if (n == 2) {
            return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
        }

        float det = 0;
        float** minor = new float* [n - 1];
        for (int i = 0; i < n - 1; ++i) {
            minor[i] = new float[n - 1];
        }

        for (int col = 0; col < n; ++col) {
            getMinor(matrix, minor, 0, col, n);
            float cofactor = (col % 2 == 0 ? 1 : -1) * matrix[0][col] * determinant(minor, n - 1);
            det += cofactor;
        }

        // Clean up
        for (int i = 0; i < n - 1; ++i) {
            delete[] minor[i];
        }
        delete[] minor;

        return det;
    }

};

// Overloaded input/output operators
ostream& operator<<(ostream& os, const Matrix& m) {
    for (int i = 0; i < m.getrow(); ++i) {
        for (int j = 0; j < m.getcol(); ++j) {
            os << m.getmat()[i][j] << " ";
        }
        os << endl;
    }
    return os;
}

istream& operator>>(istream& is, Matrix& m) {
    for (int i = 0; i < m.getrow(); ++i) {
        for (int j = 0; j < m.getcol(); ++j) {
            is >> m.getmat()[i][j];
        }
    }
    return is;
}int main() {
    // Take input for matrix dimensions
    int rows, cols;
    cout << "Enter the number of rows and columns for Matrix A: ";
    cin >> rows >> cols;
    Matrix A(rows, cols);

    cout << "Enter elements of Matrix A:" << endl;
    cin >> A;

    // Display Matrix A
    cout << "Matrix A:" << endl;
    cout << A;

    // Take input for Matrix B with same dimensions
    cout << "Enter the number of rows and columns for Matrix B: ";
    cin >> rows >> cols;
    Matrix B(rows, cols);

    cout << "Enter elements of Matrix B:" << endl;
    cin >> B;

    // Display Matrix B
    cout << "Matrix B:" << endl;
    cout << B;

    // Matrix Addition
    if (A.getrow() == B.getrow() && A.getcol() == B.getcol()) {
        Matrix C = A + B;
        cout << "Matrix A + B:" << endl;
        cout << C;

        // Matrix Subtraction
        Matrix D = A - B;
        cout << "Matrix A - B:" << endl;
        cout << D;
    }
    else {
        cout << "Matrix dimensions must be the same for addition and subtraction." << endl;
    }

    // Scalar Multiplication
    float scalar;
    cout << "Enter a scalar value to multiply Matrix A: ";
    cin >> scalar;
    Matrix E = A * scalar;
    cout << "Matrix A * " << scalar << ":" << endl;
    cout << E;

    // Matrix Multiplication

    if (A.getcol() == B.getrow()) {
        Matrix G = A * B;
        cout << "Matrix A * B:" << endl;
        cout << G;
    }
    else {
        cout << "Matrices cannot be multiplied due to incompatible dimensions." << endl;
    }

    // Transpose of Matrix A
    Matrix H = A.transpose();
    cout << "Transpose of Matrix A:" << endl;
    cout << H;

    // Inverse of Matrix A
    Matrix I = A.inverse();
    if (I.getrow() != 0 && I.getcol() != 0) {  // Check if the inverse operation was successful
        cout << "Inverse of Matrix A:" << endl;
        cout << I;
    }
    else {
        cout << "Matrix A is not invertible." << endl;
    }

    // Determinant of Matrix A
    float det = A.determinant();
    cout << "Determinant of Matrix A: " << det << endl;
    // Define matrix size
    const int size = 3; // Adjust as needed

    // Initialize a Matrix object (assuming it's a square matrix for QR decomposition and eigenvalue calculation)
    Matrix matrix(size, size);

    // Input matrix elements
    cout << "Enter elements of the matrix:" << endl;
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            cin >> matrix.getmat()[i][j];
        }
    }

    // Compute eigenvalues
    float* eigenvals = matrix.eigenvalues();
    cout << "Eigenvalues of the matrix:" << endl;
    for (int i = 0; i < size; ++i) {
        cout << eigenvals[i] << " ";
    }
    cout << endl;

    // Clean up dynamically allocated memory
    delete[] eigenvals;
    return 0;
}
