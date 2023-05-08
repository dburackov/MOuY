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


def get_base_matrix(a_matrix, b_base_plan):
    result = []
    
    for i in range(len(b_base_plan)):
        curr = []
        for j in b_base_plan:
            curr.append(a_matrix[i][j])
        result.append(curr)
    
    return result
