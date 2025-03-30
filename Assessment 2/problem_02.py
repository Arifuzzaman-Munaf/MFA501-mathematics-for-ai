# This class will be for all the operation needed for polynomial
class Polynomial:
    def __init__(self, coefficients):
        self.coefficients = coefficients
        self.degree = len(coefficients)

    def get_degree(self):
        return self.degree  # return degree of polynomial

    def get_coefficients(self):
        return self.coefficients  # return coefficients of polynomial

    def __add__(self, other):
        """
        add two polynomials
        :param other: polynomial instance
        :return: a polynomial instance after adding the coefficients
        """
        max_len = max(self.degree, other.degree)  # get the size of polynomial with more coefficients
        result = [0] * max_len  # set the size of new polynomial
        for i in range(max_len):
            result_c1 = self.coefficients[i] if i < self.degree else 0  # use 0 if index exceeds self's degree
            result_c2 = other.coefficients[i] if i < other.degree else 0
            result[i] = result_c1 + result_c2  # Sum the coefficients from both polynomials stored in result

        return Polynomial(result)

    def __sub__(self, other):
        """
        subtract the coefficients of polynomial instance 'other' from 'self'
        :param other: polynomial instance
        :return: a polynomial instance after subtracting the coefficients
        """
        other = other * Polynomial([-1])  # multiply the coefficient of 'other' with -1 to subtract
        return Polynomial((self + other).coefficients)  # we have used __add__ for subtraction as self + (-other)

    def __mul__(self, other):
        """
        Multiplies two polynomials
        :param other: polynomial instance
        :return: a polynomial instance storing multiplication between the coefficients of 'self' and 'other'
        """
        n = len(self.coefficients) + other.degree - 1  # determine the size of new polynomial
        result = [0] * n  # create a zero matrix of n

        for i in range(self.degree):
            for j in range(other.degree):
                result[i + j] += self.coefficients[i] * other.coefficients[j]

        return Polynomial(result)

    def __repr__(self):
        return f"Polynomial({self.coefficients})"  # return the list of coefficient in this format


# This class will be for all the operation needed for a nxn Matrix
class Matrix:
    def __init__(self, matrix):
        self.order = len(matrix)
        self.matrix = matrix

    def __getitem__(self, row):
        return self.matrix[row]  # return row_list of matrix denoting 'key' as row number

    def __setitem__(self, row, value):
        self.matrix[row] = value  # set the value of row

    def __mul__(self, other):
        """
        Multiplies two matrices.
        :param other: nxn matrix represented as a Matrix instance.
        :return: A new Matrix instance representing the product.
        """
        n = self.order
        # Compute the transpose of the second matrix for faster column access.
        other_T = [[other[j][i] for j in range(n)] for i in range(n)]
        # Use list comprehensions and zip to compute the dot products.
        result = [[sum(a * b for a, b in zip(self[i], other_T[j])) for j in range(n)] for i in range(n)]
        return Matrix(result)

    def get_characteristic_matrix(self):
        """
        create a polynomial matrix using A-λI formula
        :return: a polynomial matrix where coefficient of 0th and 1st position represent λ and A[i][j]
        """
        for i in range(self.order):
            for j in range(self.order):
                self[i][j] = Polynomial([-1 if i == j else 0, self[i][j]])
        return self

    def get_eigen_matrix(self, lamda):
        """
        create eigen matrix using A-λI formula
        :param lamda: a float, eigenvalue λ
        :return: a nxn matrix
        """
        n = self.order
        eigen_matrix = [[0] * n for i in range(n)]  # create zero matrix
        for i in range(n):
            for j in range(n):
                eigen_matrix[i][j] = self[i][j] - (lamda if i == j else 0)  # subtracting eigenvalue from diagonal
        return Matrix(eigen_matrix)  # return eigen matrix as a Matrix instance

    def create_copy(self):
        """
        create a copy of matrix at different reference
        :return: a matrix instance
        """
        return Matrix([row[:] for row in self.matrix])


def get_char_equation(ch_matrix):
    """
    a recursive function to generate characteristic equation
    :param ch_matrix: characteristic matrix with the coefficients of λ
    :return: a polynomial instance representing characteristic equation
    """
    order = ch_matrix.order
    if order == 2:
        # Base case for 2x2 matrices
        return (ch_matrix[0][0] * ch_matrix[1][1]) - (ch_matrix[0][1] * ch_matrix[1][0])  # perform a00*b11 - a01* b10

    poly = Polynomial([])  # initialised an empty polynomial
    sub_rows = ch_matrix[1:]  # Fetch all rows except the first row
    for i, value in enumerate(ch_matrix[0]):  # Expand through elements of the first row
        minor = Matrix(
            [row[:i] + row[i + 1:] for row in sub_rows])  # Create sub-matrix (minor) by removing the i-th column
        term = value * get_char_equation(minor)  # create equation for minor and multiply with ch_matrix[0][i]

        # calculates the sign using (-1)^row+column; row is always 1
        # if sign is negative, subtract the polynomial else add
        poly = poly - term if (-1) ** (i + 1) == -1 else poly + term
    return poly


def get_companion_matrix(polynomial, n):
    """
    creates a transposed companion matrix
    :param polynomial: a polynomial instance
    :param n: order of nxn matrix
    :return: a matrix instance, containing transposed companion matrix
    """
    coeffs = polynomial.get_coefficients()  # [a0, a1, ..., an] for a0*x^n + a1*x^(n-1) + ... + an,
    # create 1st row: [-a1, -a2, ..., -an]
    first_row = [-coeffs[i] for i in range(1, n + 1)]
    # create Remaining rows with sub-diagonal ones
    comp_matrix = [first_row]
    for i in range(1, n):
        row, row[i - 1] = [0] * n, 1  # create a list of zeros of n size and set diagonal = 1
        comp_matrix.append(row)
    return Matrix(comp_matrix)


def qr_decomposition(A):
    """
    # perform QR decomposition using classical Gram–Schmidt to make A = Q * R
    :param A: nxn companion matrix
    :return: two instances of Matrix consisting Q and R matrix after decomposition
    """
    n = A.order
    Q = [[0] * n for j in range(n)]  # Initialize Q as an n x n zero matrix
    R = [[0] * n for j in range(n)]  # Initialize R as n x n zero matrix.
    # Process each column j of A.
    for j in range(n):
        # Copy the j-th column of A into vector v.
        v = [A[i][j] for i in range(n)]

        # Subtract the projections of v on all previously computed q vectors (columns of Q).
        for i in range(j):
            # Extract the i-th column of Q as vector q.
            q = [Q[k][i] for k in range(n)]
            # Compute the dot product of q and v.
            R[i][j] = sum(q[k] * v[k] for k in range(n))
            # Subtract the projection of v onto q from v.
            for k in range(n):
                v[k] -= R[i][j] * q[k]

        # Compute the norm of the modified vector v.
        norm_v = sum(x * x for x in v) ** 0.5
        # The norm becomes the diagonal entry R[j][j].
        R[j][j] = norm_v

        if norm_v == 0:
            # If the norm is zero, then the j-th column of Q remains as zeros.
            for k in range(n):
                Q[k][j] = 0
        else:
            # Otherwise, normalize v to form the j-th column of Q.
            for k in range(n):
                Q[k][j] = v[k] / norm_v
    return Matrix(Q), Matrix(R)  # Return the matrices Q and R as Matrix instance.


def qr_algorithm(c_matrix, max_iter=1000, tol=1e-9):
    """
    Perform the QR algorithm to compute eigenvalues of a matrix.
    :param c_matrix: companion matrix of nxn order
    :param max_iter: Maximum number of iterations to perform.
    :param tol: Tolerance for convergence (based on off-diagonal elements).
    :return: List of eigenvalues
    """
    n = c_matrix.order  # Get the size of the matrix (n x n).
    for iteration in range(max_iter):  # Iterate up to the maximum number of iterations.
        Q, R = qr_decomposition(c_matrix)  # Perform QR decomposition: c_matrix = Q * R
        A_next = R * Q  # Form the next matrix in the sequence: A_{k+1} = R * Q

        # Check for convergence by calculating the sum of absolute values
        off_diag = 0
        for i in range(n):  # Iterate over rows
            for j in range(n):  # Iterate over columns
                if i != j:  # Only consider off-diagonal elements
                    off_diag += abs(A_next[i][j])

        # If the sum of off-diagonal elements is below the tolerance, we assume convergence.
        if off_diag < tol:
            c_matrix = A_next  # Update c_matrix to the converged matrix.
            break
        c_matrix = A_next  # Update c_matrix for the next iteration.
    return [c_matrix[i][i] for i in range(n)]  # list of diagonal elements of the final matrix.


# The np.roots-like function:
def get_eigenvalues(polynomial, max_iter=1000, tol=1e-9):
    """
    Find the roots of a polynomial (real or complex) with coefficients given
    in descending order, using the companion matrix and QR algorithm.
    :param polynomial: the characteristic equation as Polynomial instance
    :param max_iter: Maximum number of iterations to perform.
    :param tol: Tolerance for convergence (used in the QR algorithm).
    :return: A set of approximated eigenvalues, rounded to 2 decimal places.
    """
    n = polynomial.get_degree() - 1  # set the number of eigenvalues
    if n < 1:
        return []  # No roots for constant polynomial.
    com_matrix = get_companion_matrix(polynomial, n)  # get companion matrix
    eigenvalues = qr_algorithm(com_matrix, max_iter, tol)  # perform QR algorithm on companion matrix
    return set(round(eigenvalue, 2) for eigenvalue in eigenvalues)


def gaussian_elimination(mat, tol=1e-8):
    """
    Perform Gaussian elimination on a nxn matrix.
    :param mat: nxn matrix
    :param tol: Tolerance for considering a value as zero.
    :return: The matrix in row echelon form.
    """
    m = mat.create_copy()  # create a copy of the matrix object.
    rows, cols = m.order, len(m[0])  # Get the number of rows and columns in the matrix.
    pivot_row = 0  # Keep track of the current pivot row.
    for col in range(cols - 1):
        pivot = None  # Initialize pivot as None
        for r in range(pivot_row, rows):  # Start searching from the current pivot row.
            if abs(m[r][col]) > tol:  # Check if the element is larger than the tolerance.
                pivot = r  # Set pivot to the row index where a valid pivot is found.
                break  # Stop searching once a valid pivot is found.

        if pivot is None:  # If no valid pivot is found, move to the next column.
            continue

        m[pivot_row], m[pivot] = m[pivot], m[
            pivot_row]  # Swap the current pivot row with the row containing the pivot element.
        factor = m[pivot_row][col]  # Normalize the pivot row by dividing all elements by the pivot value.
        m[pivot_row] = [val / factor for val in m[pivot_row]]  # Normalize row.

        # Eliminate entries below the pivot by subtracting multiples of the pivot row.
        for r in range(pivot_row + 1, rows):  # Iterate over rows below the pivot row.
            factor = m[r][col]  # Factor to eliminate current column entry.
            m[r] = [m[r][c] - factor * m[pivot_row][c] for c in range(cols)]  # Row operation.

        pivot_row += 1  # Move to the next pivot row for subsequent columns.
        if pivot_row == rows:  # If all rows have been processed, stop early.
            break

    return m


def get_nullspace(eigen_matrix, tol=1e-8):
    """
    :param eigen_matrix: a matrix instance of nxn matrix
    :param tol: Tolerance for considering a value as zero.
    :return: list of lists as nullspace
    """
    n = eigen_matrix.order  # Get the dimension (order) of the matrix.
    # Create an augmented matrix by appending a zero to each row
    aug_matrix = Matrix([row[:] + [0] for row in eigen_matrix.matrix])
    rref = gaussian_elimination(aug_matrix)  # get the row-reduced echelon form (RREF).
    pivot_cols = []  # Initialize list to store indices of pivot columns.

    for row in rref:  # Iterate over each row in the RREF matrix.
        for j in range(n):  # Loop over columns.
            if abs(row[j]) > tol:  # If the element is significantly nonzero, it's a pivot.
                pivot_cols.append(j)  # store the pivot column index.
                break  # Move to the next row
    free_vars = [j for j in range(n) if j not in pivot_cols]  # Identify free variable indices

    if not free_vars:  # If there are no free variables,
        return [[0] * n]  # Return the trivial solution.

    basis = []  # Initialize list to hold basis vectors.
    for free in free_vars:
        sol = [0] * n  # Start with a solution vector of zeros.
        sol[free] = 1  # Set the free variable to 1 for this basis vector.
        for row in rref[::-1]:  # Perform back substitution: iterate over rows in reverse.
            # Find the pivot column in the current row (first nonzero element).
            pivot_col = next((j for j in range(n) if abs(row[j]) > tol), None)
            if pivot_col is None:  # If no pivot is found in this row, continue.
                continue
            # Compute the sum of the contributions of the known variables in the row.
            s = sum(row[j] * sol[j] for j in range(pivot_col + 1, n))
            sol[pivot_col] = -s  # Solve for the pivot variable.
        basis.append(sol)  # Add the computed solution vector to the basis list.
    return basis  # Return the complete nullspace basis.


def get_eigenspace(matrix, eigenvalue):
    """
    :param matrix: instance nxn matrix
    :param eigenvalue: integer representing eigenvalue
    :return: list of list (list of vectors)
    """
    eigen_matrix = matrix.get_eigen_matrix(eigenvalue)
    return get_nullspace(eigen_matrix)


matrix = Matrix([[5, 4, 2],
                 [4, 5, 2],
                 [2, 2, 2],
                 ])
print(f"Actual Matrix: {matrix.matrix}\n")

mat = matrix.create_copy()  # create a copy of instance of matrix
characteristic_matrix = mat.get_characteristic_matrix()

characteristic_polynomial = get_char_equation(characteristic_matrix)
print(f"Characteristic Polynomial = {characteristic_polynomial}\n")

eigenvalues = get_eigenvalues(characteristic_polynomial)
print(f"Eigenvalues of {characteristic_polynomial} = {eigenvalues}\n")

for lamb_da in eigenvalues:
    basis = get_eigenspace(matrix, lamb_da)
    print(f"Eigen-space basis for λ = {lamb_da} : ")
    for vec in basis:
        print([round(x, 4) for x in vec])  # Format the vector to 4 decimal point.
    print()
