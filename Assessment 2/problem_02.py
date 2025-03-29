class Polynomial:
    def __init__(self, coefficients):
        self.coefficients = coefficients
        self.order = len(coefficients)

    def get_order(self):
        return self.order

    def get_coefficients(self):
        return self.coefficients

    def __add__(self, other):
        max_len = max(self.order, len(other.coefficients))
        result = [0] * max_len

        for i in range(max_len):
            coeff1 = self.coefficients[i] if i < self.order else 0
            coeff2 = other.coefficients[i] if i < len(other.coefficients) else 0
            result[i] = coeff1 + coeff2

        return Polynomial(result)

    def __sub__(self, other):
        max_len = max(self.order, len(other.coefficients))
        result = [0] * max_len

        for i in range(max_len):
            coeff1 = self.coefficients[i] if i < len(self.coefficients) else 0
            coeff2 = other.coefficients[i] if i < len(other.coefficients) else 0
            result[i] = coeff1 - coeff2

        return Polynomial(result)

    def __mul__(self, other):
        result_size = len(self.coefficients) + len(other.coefficients) - 1
        result = [0] * result_size

        for i in range(len(self.coefficients)):
            for j in range(len(other.coefficients)):
                result[i + j] += self.coefficients[i] * other.coefficients[j]

        return Polynomial(result)

    def __repr__(self):
        return f"Polynomial({self.coefficients})"





def determinant(matrix):
    """
    Compute the determinant of a matrix recursively as a polynomial.
    """
    order = len(matrix)

    if order == 2:
        # Base case for 2x2 matrices
        a0 = matrix[0][0]
        a1 = matrix[0][1]
        b0 = matrix[1][0]
        b1 = matrix[1][1]

        return (a0 * b1) - (a1 * b0)

    poly = Polynomial([0])
    sub_rows = matrix[1:]  # Fetch all rows except the first row
    for i, value in enumerate(matrix[0]):  # Expand through elements of the first row
        minor = [row[:i] + row[i + 1:] for row in sub_rows]  # Create sub-matrix (minor) by removing the i-th column

        # Compute sign using (-1)^(row+col), here row=0 so only column matters
        term = value * determinant(minor)

        if (-1) ** (i + 1) == -1:
            poly -= term
        else:
            poly += term

    return poly

def convAI(matrix):
    """
    Convert a matrix into its characteristic polynomial form.
    """
    order = len(matrix)
    for i in range(order):
        for j in range(order):
            if i == j:
                matrix[i][i] = Polynomial([-1, matrix[i][i]])
            else:
                matrix[i][j] = Polynomial([0, matrix[i][j]])

    return matrix
# Example usage


def companion_matrix(coeffs):
    """
    Given polynomial coefficients in descending order:
         [a0, a1, ..., an] for a0*x^n + a1*x^(n-1) + ... + an,
    construct the companion matrix (of size n x n).
    """
    n = len(coeffs) - 1
    a0 = coeffs[0]
    # First row: [-a1/a0, -a2/a0, ..., -an/a0]
    first_row = [-coeffs[i] / a0 for i in range(1, len(coeffs))]
    # Remaining rows: subdiagonal ones
    C = [first_row]
    for i in range(1, n):
        row = [0] * n
        row[i - 1] = 1
        C.append(row)
    return C


# Matrix multiplication (square matrices)
def matrix_multiply(A, B):
    n = len(A)
    result = [[0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            s = 0
            for k in range(n):
                s += A[i][k] * B[k][j]
            result[i][j] = s
    return result


# QR decomposition using classical Gram–Schmidt
def qr_decomposition(A):
    n = len(A)
    # Initialize Q as an n x n zero matrix, and R as n x n zero matrix.
    Q = [[0 for _ in range(n)] for _ in range(n)]
    R = [[0 for _ in range(n)] for _ in range(n)]
    # Process column by column.
    for j in range(n):
        # Copy j-th column of A into vector v.
        v = [A[i][j] for i in range(n)]
        # Subtract projections on previous q vectors.
        for i in range(j):
            # Compute dot product of q_i (column i of Q) with v.
            q = [Q[k][i] for k in range(n)]
            R[i][j] = sum(q[k] * v[k] for k in range(n))
            for k in range(n):
                v[k] -= R[i][j] * q[k]
        # Compute norm of v.
        norm_v = sum(x * x for x in v) ** 0.5
        R[j][j] = norm_v
        if norm_v == 0:
            # If v is zero, leave Q's column as zeros.
            for k in range(n):
                Q[k][j] = 0
        else:
            for k in range(n):
                Q[k][j] = v[k] / norm_v
    return Q, R


# The basic QR algorithm (without shifts) for eigenvalue computation.
def qr_algorithm(A, max_iter=1000, tol=1e-9):
    n = len(A)
    # Make a copy of A.
    Ak = [row[:] for row in A]
    for iteration in range(max_iter):
        Q, R = qr_decomposition(Ak)
        # Form next matrix: A_{k+1} = R * Q
        Ak_next = matrix_multiply(R, Q)
        # Check convergence: sum of off-diagonal absolute values.
        off_diag = 0
        for i in range(n):
            for j in range(n):
                if i != j:
                    off_diag += abs(Ak_next[i][j])
        if off_diag < tol:
            Ak = Ak_next
            break
        Ak = Ak_next
    # The eigenvalues are the diagonal elements.
    return [Ak[i][i] for i in range(n)]


# The np.roots-like function:
def np_roots(polynomial, max_iter=1000, tol=1e-9):
    """
    Find the roots of a polynomial (real or complex) with coefficients given
    in descending order, using the companion matrix and QR algorithm.

    Returns:
         list of approximated eigenvalues (roots).
    """
    n = polynomial.get_order() - 1
    if n < 1:
        return []  # No roots for constant polynomial.
    C = companion_matrix(polynomial.coefficients)
    eigenvalues = qr_algorithm(C, max_iter, tol)
    return set(round(eigenvalue, 2) for eigenvalue in eigenvalues)


matrix = [
    [5, 4, 2,4],
    [4, 5, 2,3],
    [2, 2, 2,1],
    [4, 3, 1,2]
    ]

# Convert the matrix to its characteristic polynomial form
characteristic_matrix = convAI(matrix)
print(f"Matrix after conversion: {characteristic_matrix}")

# Compute the characteristic polynomial
characteristic_polynomial = determinant(characteristic_matrix)
print(f"Characteristic Polynomial: {characteristic_polynomial.coefficients}")
roots = np_roots(characteristic_polynomial)

print(roots)

