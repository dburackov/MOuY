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
        for i in range(len(b_base_plan)):
            b_base_plan[i] -= 1
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
        return x_base_plan, b_base_plan
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


def simplex_method_init(c_vec, a_mtx, b_vec):
    n = len(a_mtx[0])
    m = len(a_mtx)
    #1
    for i in range(m):
        if b_vec[i] < 0:
            b *= -1
            a_mtx[i] = (np.array(a_mtx[i]) * -1).tolist()
    #2
    c_2 = [0] * n + [-1] * m

    a_2 = [a_mtx[i].copy() for i in range(m)]
    id_mtx = np.eye(m)
    for i in range(m):
        for j in range(m):
            a_2[i].append(id_mtx[i][j])
    #3
    x_init_plan = [0] * n + b_vec
    b_init_plan = [(n + i + 1) for i in range(m)]
    #4
    x_init_plan, b_init_plan = simplex_method_main(c_2, a_2, x_init_plan, b_init_plan)
    #5
    for i in range(m):
        if x_init_plan[n + i] != 0:
            raise Exception("задача не совместима")
    #6
    x_base_plan = x_init_plan[:n]
    #7
    while max(b_init_plan) >= n:
        #8
        j_k = max(b_init_plan)
        k = b_init_plan.index(j_k)
        #7
        a_2_inv = np.linalg.inv(get_base_matrix(a_2, b_init_plan)) 

        flag = False
        for j in range(n):
            if j not in b_init_plan:
                l = a_2_inv.dot(np.array(a_2)[:, j])
                if l[k] != 0:
                    b_init_plan[k] = j
                    flag = True
                    break

        if not flag:  
            i = j_k - n
            a_2.pop(i)
            a_mtx.pop(i)
            b_vec.pop(i)
            b_init_plan.remove(j_k)


    return x_base_plan, b_init_plan, a_mtx, b_vec


def lab3():
    a = [[1, 1, 1], [2, 2, 2]]
    c = [1, 0, 0]
    b = [0, 0]

    print(simplex_method_init(c, a, b))

    a = [[0, -1, 1, 1, 0],
                [-5, 1, 1, 0, 0],
                [-8, 1, 2, 0, -1]]
    c = [-3, 1, 4, 0, 0]
    b = [1, 2, 3]

    print(simplex_method_init(c, a, b))


def lab2():
    c = [1, 1, 0, 0, 0]
    x = [0, 0, 1, 3, 2]
    a = [[-1, 1, 1, 0, 0], [1, 0, 0, 1, 0], [0, 1, 0, 0, 1]]
    b = [3, 4, 5]
    print(simplex_method_main(c, a, x, b))
    print(x)
    print(b)


if __name__ == '__main__':
    lab3()
    



