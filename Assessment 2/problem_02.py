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


class Matrix:
    def __init__(self, matrix):
        self.order = len(matrix)
        self.matrix = matrix

    def __getitem__(self, row):
        return self.matrix[row]

    def __setitem__(self, row, value):
        self.matrix[row] = value

    def __mul__(self, other):
        n = self.order
        result = [[0] * n for i in range(n)]
        for i in range(n):
            for j in range(n):
                s = 0
                for k in range(n):
                    s += self[i][k] * other[k][j]
                result[i][j] = s
        return Matrix(result)

    def get_characteristic_matrix(self):
        for i in range(self.order):
            for j in range(self.order):
                self[i][j] = Polynomial([-1 if i == j else 0, self[i][j]])
        return self

    def get_eigen_matrix(self, lamda):
        n = self.order
        eigen_matrix = [[0] * n for i in range(n)]
        for i in range(n):
            for j in range(n):
                eigen_matrix[i][j] = self[i][j] - (lamda if i == j else 0)
        return Matrix(eigen_matrix)


def get_polynomial(ch_matrix):
    """
    Compute the determinant of a matrix recursively as a polynomial.
    """
    order = ch_matrix.order

    if order == 2:
        # Base case for 2x2 matrices
        return (ch_matrix[0][0] * ch_matrix[1][1]) - (ch_matrix[0][1] * ch_matrix[1][0])

    poly = Polynomial([0])
    sub_rows = ch_matrix[1:]  # Fetch all rows except the first row
    for i, value in enumerate(ch_matrix[0]):  # Expand through elements of the first row
        minor = Matrix(
            [row[:i] + row[i + 1:] for row in sub_rows])  # Create sub-matrix (minor) by removing the i-th column

        # Compute sign using (-1)^(row+col), here row=0 so only column matters
        term = value * get_polynomial(minor)

        poly = poly - term if (-1) ** (i + 1) == -1 else poly + term
    return poly


def get_companion_matrix(polynomial, n):
    """
    Given polynomial coefficients in descending order:
         [a0, a1, ..., an] for a0*x^n + a1*x^(n-1) + ... + an,
    construct the companion matrix (of size n x n).
    """
    coeffs = polynomial.get_coefficients()
    a0 = coeffs[0]
    # First row: [-a1/a0, -a2/a0, ..., -an/a0]
    first_row = [-coeffs[i] / a0 for i in range(1, n + 1)]
    # Remaining rows: subdiagonal ones
    comp_matrix = [first_row]
    for i in range(1, n):
        row = [0] * n
        row[i - 1] = 1
        comp_matrix.append(row)
    return Matrix(comp_matrix)


# QR decomposition using classical Gram–Schmidt
def qr_decomposition(matrix):
    n = matrix.order
    # Initialize Q as an n x n zero matrix, and R as n x n zero matrix.
    Q = [[0] * n for j in range(n)]
    R = [[0] * n for j in range(n)]
    # Process column by column.
    for j in range(n):
        # Copy j-th column of A into vector v.
        v = [matrix[i][j] for i in range(n)]
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
    return Matrix(Q), Matrix(R)


# The basic QR algorithm (without shifts) for eigenvalue computation.
def qr_algorithm(c_matrix, max_iter=1000, tol=1e-9):
    n = c_matrix.order
    # Make a copy of A.

    for iteration in range(max_iter):
        Q, R = qr_decomposition(c_matrix)
        # Form next matrix: A_{k+1} = R * Q
        temp_next = R * Q
        # Check convergence: sum of off-diagonal absolute values.
        off_diag = 0
        for i in range(n):
            for j in range(n):
                if i != j:
                    off_diag += abs(temp_next[i][j])
        if off_diag < tol:
            c_matrix = temp_next
            break
        c_matrix = temp_next
    # The eigenvalues are the diagonal elements.
    return [c_matrix[i][i] for i in range(n)]


# The np.roots-like function:
def get_eigenvalues(polynomial, max_iter=1000, tol=1e-9):
    """
    Find the roots of a polynomial (real or complex) with coefficients given
    in descending order, using the companion matrix and QR algorithm.

    Returns:
         list of approximated eigenvalues (roots).
    """
    n = polynomial.get_order() - 1
    if n < 1:
        return []  # No roots for constant polynomial.
    com_matrix = get_companion_matrix(polynomial, n)
    eigenvalues = qr_algorithm(com_matrix, max_iter, tol)
    return set(round(eigenvalue, 2) for eigenvalue in eigenvalues)


def gaussian_elimination(mat, tol=1e-8):
    """
    Perform Gaussian elimination on the augmented matrix (for homogeneous system).
    Returns the matrix in row-echelon form.
    """
    m = [row[:] for row in mat]  # shallow copy each row
    rows, cols = len(m), len(m[0])
    pivot_row = 0
    for col in range(cols - 1):  # exclude augmented column
        # Find pivot in the current column
        pivot = None
        for r in range(pivot_row, rows):
            if abs(m[r][col]) > tol:
                pivot = r
                break
        if pivot is None:
            continue
        # Swap and normalize pivot row
        m[pivot_row], m[pivot] = m[pivot], m[pivot_row]
        factor = m[pivot_row][col]
        m[pivot_row] = [val / factor for val in m[pivot_row]]
        # Eliminate entries below pivot
        for r in range(pivot_row + 1, rows):
            factor = m[r][col]
            m[r] = [m[r][c] - factor * m[pivot_row][c] for c in range(cols)]
        pivot_row += 1
        if pivot_row == rows:
            break
    return m


def compute_nullspace(eigen_matrix, tol=1e-8):
    """
    Given a square matrix M (without the augmented column), compute a basis for its nullspace.
    """
    n = eigen_matrix.order
    # Append augmented zero column
    augmented_matrix = [row[:] + [0] for row in eigen_matrix.matrix]
    rref = gaussian_elimination(augmented_matrix)

    # Identify pivot columns
    pivot_cols = []
    for row in rref:
        for j in range(n):
            if abs(row[j]) > tol:
                pivot_cols.append(j)
                break
    free_vars = [j for j in range(n) if j not in pivot_cols]

    # If no free variable, return the trivial solution
    if not free_vars:
        return [[0] * n]

    basis = []
    for free in free_vars:
        sol = [0] * n
        sol[free] = 1
        # Back substitution from bottom to top
        for row in reversed(rref):
            # Find the first nonzero (pivot) entry in the row
            pivot_col = next((j for j in range(n) if abs(row[j]) > tol), None)
            if pivot_col is None:
                continue
            # Adjust pivot variable based on free variables
            s = sum(row[j] * sol[j] for j in range(pivot_col + 1, n))
            sol[pivot_col] = -s
        basis.append(sol)
    return basis


def get_eigenspace(matrix, igenvalue):
    """Compute a basis for the eigenspace corresponding to eigenvalue lam."""
    eigen_matrix = matrix.get_eigen_matrix(igenvalue)
    return compute_nullspace(eigen_matrix)


matrix = [[5, 4, 2],
          [4, 5, 2],
          [2, 2, 2],
          ]
print(f"Actual Matrix: {matrix}\n")

mat = Matrix([row[:] for row in matrix])
characteristic_matrix = mat.get_characteristic_matrix()
print(f"Characteristic Matrix = {characteristic_matrix.matrix}\n")

characteristic_polynomial = get_polynomial(characteristic_matrix)
print(f"Characteristic Polynomial = {characteristic_polynomial}\n")

eigenvalues = get_eigenvalues(characteristic_polynomial)
print(f"Eigenvalues of {characteristic_polynomial} = {eigenvalues}\n")

matrix = Matrix(matrix)
for eigenvalue in eigenvalues:
    basis = get_eigenspace(matrix, eigenvalue)
    print("Eigenvalue:", eigenvalue)
    print("Eigenspace basis:")
    for vec in basis:
        # Format the vector to a few decimals.
        formatted = [round(x, 4) for x in vec]
        print(formatted)
    print()
