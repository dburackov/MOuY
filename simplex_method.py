import numpy as np
import mtrinv


def get_base_matrix(a_matrix, b_base_plan):
    result = []
    
    for i in range(len(b_base_plan)):
        curr = []
        for j in b_base_plan:
            curr.append(a_matrix[i][j])
        result.append(curr)
    
    return result

def simplex_method_main(c_vector, a_matrix, x_base_plan, b_base_plan, iter = 1, a_base_matrix_inv = []):
    if iter == 1:
        a_base_matrix_inv = np.linalg.inv(get_base_matrix(a_matrix, b_base_plan)) 
    c_base_vector = []
    for i in b_base_plan:
        c_base_vector.append(c_vector[i])
    #print(c_base_vector, 'cB')
    u_vector = (np.matrix(c_base_vector) * a_base_matrix_inv).A1
    #print(u_vector, 'potential vector')
    marks_vector = (np.matrix(u_vector) * a_matrix).A1 - c_vector
    #print(marks_vector, 'marks vector')
    if marks_vector.dot(marks_vector < 0) == 0:
        return x_base_plan
    j0 = (marks_vector < 0).tolist().index(True) 
    #print(j0, 'negative component in marks vector')
    a_j0 = []
    for i in range(len(a_matrix)):
        a_j0.append(a_matrix[i][j0])
    z_vector = a_base_matrix_inv.dot(a_j0)
    #print(z_vector, 'z vector')
    q_vector = []
    for i in range(len(z_vector)):
        if z_vector[i] > 0:
            q_vector.append(x_base_plan[b_base_plan[i]] / z_vector[i])
        else:
            q_vector.append(np.Infinity)
    #print(q_vector, 'q vector')
    q0 = min(q_vector)
    if q0 == np.Infinity: 
        raise Exception('целевой функционал задачи не ограничен сверху на множестве допустимых планов')
    k = q_vector.index(q0)
    j_k = b_base_plan[k]
    b_base_plan[k] = j0
    x_base_plan[j0] = q0

    for i in range(len(b_base_plan)):
        if i != k:
            x_base_plan[b_base_plan[i]] = (x_base_plan[b_base_plan[i]] - q0 * z_vector[i])

    x_base_plan[j_k] = 0
    #print(x_base_plan)
    #print(b_base_plan)

    tmp = mtrinv.get_inverse_matrix(a_base_matrix_inv, np.array(a_matrix)[:, j0], k)

    return simplex_method_main(c_vector, a_matrix, x_base_plan, b_base_plan, iter + 1, np.array(tmp))


if __name__ == '__main__':
    c = [1, 1, 0, 0, 0]
    x = [0, 0, 1, 3, 2]
    a = [[-1, 1, 1, 0, 0], [1, 0, 0, 1, 0], [0, 1, 0, 0, 1]]
    b = [2, 3, 4]
    print(simplex_method_main(c, a, x, b))




