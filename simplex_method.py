import numpy as np


def get_base_matrix(a_matrix, b_base_plan):
    result = []
    
    for i in range(len(b_base_plan)):
        curr = []
        for j in b_base_plan:
            curr.append(a_matrix[i][j])
        result.append(curr)
    
    return result

def simplex_method_main(c_vector, a_matrix, x_base_plan, b_base_plan, iter = 1):
    #1 step
    if iter == 1:
        pass
    else:
        pass
    a_base_matrix = get_base_matrix(a_matrix, b_base_plan)
    a_base_matrix_inv = np.linalg.inv(a_base_matrix) 

    #2 step
    c_base_vector = []
    for i in b_base_plan:
        c_base_vector.append(c_vector[i])
    #print(c_base_vector, 'cB')
    #3 step
    u_vector = (np.matrix(c_base_vector) * a_base_matrix_inv).A1
    #print(u_vector, 'potential vector')
    #4 step
    #marks_vector = a_matrix.dot(u_vector) - c_vector
    #marks_vector = np.squeeze(np.asarray(np.matrix(u_vector) * a_matrix)) - c_vector
    marks_vector = (np.matrix(u_vector) * a_matrix).A1 - c_vector
    #print(marks_vector, 'marks vector')
    #5 step
    if marks_vector.dot(marks_vector < 0) == 0:
        return x_base_plan
    #6 step
    j0 = (marks_vector < 0).tolist().index(True) 
    #print(j0, 'negative component in marks vector')
    #7 step
    a_j0 = []
    for i in range(len(a_matrix)):
        a_j0.append(a_matrix[i][j0])
    z_vector = a_base_matrix_inv.dot(a_j0)
    #print(z_vector, 'z vector')
    #8 step
    q_vector = []
    for i in range(len(z_vector)):
        if z_vector[i] > 0:
            q_vector.append(x_base_plan[b_base_plan[i]] / z_vector[i])
        else:
            q_vector.append(np.Infinity)
    #print(q_vector, 'q vector')
    #9 step
    q0 = min(q_vector)
    #10 step
    if q0 == np.Infinity: 
        Exception('целевой функционал задачи не ограничен сверху на множестве допустимых планов')
    #11 step
    k = q_vector.index(q0)
    j_k = b_base_plan[k]
    #12 step
    b_base_plan[k] = j0
    #13 step
    x_base_plan[j0] = q0

    for i in range(len(b_base_plan)):
        if i != k:
            x_base_plan[b_base_plan[i]] = (x_base_plan[b_base_plan[i]] - q0 * z_vector[i])

    x_base_plan[j_k] = 0
    #print(x_base_plan)
    #print(b_base_plan)
    return simplex_method_main(c_vector, a_matrix, x_base_plan, b_base_plan, iter + 1)


if __name__ == '__main__':
    c = [1, 1, 0, 0, 0]
    x = [0, 0, 1, 3, 2]
    a = [[-1, 1, 1, 0, 0], [1, 0, 0, 1, 0], [0, 1, 0, 0, 1]]
    b = [2, 3, 4]
    print(simplex_method_main(c, a, x, b))




