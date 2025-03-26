sign = 1     # determine the sign of value, changes sign each time any row is swapped.


def swap_row(matrix, candidate_row=1):
    """
    :param matrix: list of lists representing the nxn matrix
    :param candidate_row: rows to swap with 1st row, starting from 1
    :return list: list of lists after swapping with nth row having non-zero pivot.
    """
    if candidate_row >= len(matrix):
        return matrix  # No candidate row found; matrix remains unchanged.
    if matrix[candidate_row][0] != 0:
        # Swap rows 0 and candidate.
        matrix[0], matrix[candidate_row] = matrix[candidate_row], matrix[0]  # swapping rows with 0th row
        return matrix

    return swap_row(matrix, candidate_row + 1)  # continue recursive search until gets a non-zero pivot in 0th row


def gaussian_elimination(matrix, col_index=1):
    """
    :param matrix: list of lists representing the nxn matrix
    :param col_index: index of column to perform elimination, starting from 1
    :return list: list of list after elimination.
    """
    if col_index >= len(matrix[0]):  # base case 1 : elimination performed on all columns
        return matrix

    # base case 2
    # if the 0th row pivot is zero and only zero even after swapping, then no elimination needed
    if matrix[0][0] == 0:
        matrix = swap_row(matrix, 1)
        if matrix[0][0] == 0:
            return matrix  # terminate and return the actual matrix if 0th zero pivot
        else:
            global sign
            sign = -sign  # got non-zero pivot after swap, hence changed the sign for determinant

    # Calculate the reduction factor to perform Invariance
    reduction_factor = matrix[0][col_index] / matrix[0][0]

    # For every row, eliminate the col_index-th element in the first row.
    # This operation will zero out all matrix[0][col_index].
    for i in range(len(matrix)):
        matrix[i][col_index] -= reduction_factor * matrix[i][0]  # perform Invariance: b-ka
        print(matrix)

    return gaussian_elimination(matrix, col_index + 1)  # continue recursion until any base case executes


def determinant(matrix):
    """
    :param matrix: n x n matrix as a list of lists.
    :return float: The determinant of the given n x n matrix.
    """
    order = len(matrix)  # get the order of matrix

    if order == 1:  # base case 1: 1 x 1 square matrix
        return matrix[0][0]  # return the only value in 1x1 matrix

    if order == 2:  # base case 2: 2 x 2 square matrix
        return sign * (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0])  # terminates the recursion

    reduced_matrix = gaussian_elimination(matrix)  # perform gaussian elimination
    sub_rows = reduced_matrix[1:]  # fetch all rows except first row as this row will be used for expanding
    minor = [row[1:] for row in sub_rows]  # creates sub-matrix(minor) by removing i-th column

    return reduced_matrix[0][0] * determinant(minor)


def generate_matrix(order):
    """
    :param order: the order of nxn matrix
    :return list: list of lists representing the matrix
    """
    random_elements = [i for i in range(1, order + 1)]  # create a list starting from 1 - order
    mat = []  # empty list to store matrix
    for i in range(order):
        mat.append(random_elements)
        # random formula to create random values between 1 and 49.
        random_elements = [int(i ** 2.71) % 50 for i in random_elements]
    return mat


n = int(input("What is the order of the matrix? : "))
matrix = generate_matrix(n)  # creates a square matrix of order n

# matrix = [[2, 2, 4, 5], [3, 1, 9, 2], [5, 2, 12, 2], [11, 7, 7, 9]]
print(f"Matrix = {matrix}")  # printing the matrix of order n

print(f'The determinant = {determinant(matrix)}')
