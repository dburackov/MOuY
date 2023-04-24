import numpy as np


def get_inverse_matrix(inv_matrix, x_vector, index):
    n = len(inv_matrix)

    l_vector = np.dot(inv_matrix, x_vector)

    if l_vector[index] == 0:
        raise Exception('Matrix is irreversible')
    
    l1_vector = l_vector.copy()
    l1_vector[index] = -1
    l2_vector = np.dot(l1_vector, -1 / l_vector[index])

    q_matrix = np.eye(n)
    for j in range(n):
        q_matrix[j][index] = l2_vector[j]
    
    result = []
    for i in range(n):
        curr = []
        for j in range(n): 
            if i != index:
                curr.append(inv_matrix[i][j] * q_matrix[i][i] + inv_matrix[index][j] * q_matrix[i][index])
            else:
                curr.append(inv_matrix[index][j] * q_matrix[i][index])
        result.append(curr)

    return result 


def test():
    matrix_example = [[1, 0, 5],
                      [2, 1, 6],
                      [3, 4, 0]]
    inverse_matrix_example = [[-24, 20, -5],
                              [18, -15, 4],
                              [5, -4, 1]]
    x_vector_example = [2, 2, 2]
    i = 2
    print(get_inverse_matrix(
                                   inverse_matrix_example,
                                   x_vector_example, i - 1))
    matrix_example2 = [[2, 5, 7],
                       [6, 3, 4],
                       [5, -2, -3]]
    inverse_matrix_example2 = [[1, -1, 1],
                               [-38, 41, -34],
                               [27, -29, 24]]
    x_vector_example2 = [2, 1, 3]
    i2 = 1
    print(get_inverse_matrix(
                                   inverse_matrix_example2,
                                   x_vector_example2, i2 - 1))

    matrix_example3 = [[1,-1,3,2,4],
                       [2,-3,2,-2,-6],
                       [3,-5,1,4,2],
                       [4,-2,-1,-3,-1],
                       [5,-1,0,1,2]]
    inverse_matrix_example3 = [[-21/650, 4/75, -113/1950, -47/975, 84/325],
                                [-11/325, 1/25, -61/325, -68/325, 88/325],
                                [161/650, 11/75, -217/1950, -73/975, 6/325],
                                [-97/650, 1/25, 53/650, -93/325, 63/325],
                                [9/65, -2/15, 2/195, 31/195, -7/65]]
    x_vector_example3 = [2, 1, 3,-1,2]
    i3 = 4
    print(np.matrix(get_inverse_matrix(
                                   inverse_matrix_example3,
                                   x_vector_example3, i3 - 1)))

if __name__ == '__main__':
    test()