"""
The function takes n X n matrix and fetches the sub-matrices using recursion.
the recursion continues until the base case of 2x2 is reached.
"""
def determinant(matrix):
    """
    :param matrix: n x n matrix as a list of lists.
    :return float: The determinant of the given n x n matrix.
    """
    order = len(matrix)     # get the order of matrix

    if order == 1:    # base case for 1 x 1 square matrix
        return matrix[0][0]   # return the only value in 1x1 matrix

    if order == 2:   # base case for 2 x 2 square matrix
        return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]  # terminates the recursion

    result = 0              # initialized for storing determinant value
    sub_rows = matrix[1:]   # fetch all rows except first row as this row will be used for expanding
    for i, value in enumerate(matrix[0]):                   # expands through elements of 1st row
        minor = [row[:i] + row[i + 1:] for row in sub_rows]  # creates sub-matrix(minor) by removing i-th column

        # first, compute the sign using (-1)^(col+row) formula,here row = 0,thus only column is used.
        # secondly, current entity is multiplied with sign
        # finally, operate recursion to get the determinant of minor.
        result += (-1)**i * value * determinant(minor)
    return result


matrix = [[1, 2, 3],
          [0, 4, 5],
          [1, 0, 6]]
det = determinant(matrix)

print(f'The determinant = {determinant(matrix)}')







