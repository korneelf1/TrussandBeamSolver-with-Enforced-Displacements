from audioop import cross
import unittest
from MeshGenerator import *
from ExternalForcesMomentsConstraints import Force, Moment, Support
import numpy.testing as nptest


class Test_Rod(unittest.TestCase):
    def test_Rod(self):
        rod = Rod((0,0),(1,0),cross_area_initial=100e-6)
        self.assertEqual(rod.A_i,100e-6)
        rod.A_i = 1
        self.assertEqual(rod.A_i,1)


class Test_Mesh(unittest.TestCase):
    def test_search_node(self):
        mesh1 = Mesh((0,0),(0,10),(10,0),(10,10))
        node_coords = [(0,0),(1,0),(0,1),(1,1)]
        node_connections = [(0,1),(1,2),(0,2),(1,3),(2,3)]
        
        mesh1.manual_mesh(node_coords, node_connections, seed = 1)

        self.assertEqual(mesh1.search_node((0,0)).coord, (0,0))
        self.assertIsNone(mesh1.search_node((10,0)))


    def test_add_external_force(self):
        mesh2 = Mesh((0,0),(0,10),(10,0),(10,10))
        node_coords = [(0,0),(1,0),(0,1),(1,1)]
        node_connections = [(0,1),(1,2),(0,2),(1,3),(2,3)]
        
        mesh2.manual_mesh(node_coords, node_connections, seed = 1)

        force = Force(100,0,(1,1))
        mesh2.add_external_force(force)
        self.assertEqual(len(mesh2.forces_lst),1)
        self.assertEqual(mesh2.forces_lst[0].magn, 100)
        self.assertEqual(mesh2.forces_lst[0].angle, 0)
        self.assertEqual(mesh2.forces_lst[0].node, (1,1))
       
    def test_add_external_moment(self):
        mesh3 = Mesh((0,0),(0,10),(10,0),(10,10))
        node_coords = [(0,0),(1,0),(0,1),(1,1)]
        node_connections = [(0,1),(1,2),(0,2),(1,3),(2,3)]
        
        mesh3.manual_mesh(node_coords, node_connections, seed = 1)

        moment = Moment(100,(1,1))
        mesh3.add_external_moment(moment)
        self.assertEqual(len(mesh3.moments_lst),1)
        self.assertEqual(mesh3.moments_lst[0].magn, 100)
        self.assertEqual(mesh3.moments_lst[0].node, (1,1))
       
    def test_add_support(self):
        mesh4 = Mesh((0,0),(0,10),(10,0),(10,10))
        node_coords = [(0,0),(1,0),(0,1),(1,1)]
        node_connections = [(0,1),(1,2),(0,2),(1,3),(2,3)]
        
        mesh4.manual_mesh(node_coords, node_connections, seed = 1)

        support = Support((0,0),'f',orientation=0)
        mesh4.add_support(support)
        self.assertEqual(len(mesh4.support_lst),1)
        self.assertEqual(mesh4.support_lst[0].type, 'f')
        self.assertEqual(mesh4.support_lst[0].node, (0,0))
       
        support = Support((0,0),'d',orientation=0)
        self.assertWarns(UserWarning)
        mesh4.add_support(support)
        self.assertWarns(UserWarning)

        support = Support((0,10),'f',orientation=0)
        self.assertWarns(UserWarning)

    def test_element_matrix_rod(self):
        E = 50e7
        A = .1
        start = Node(0,0)
        end = Node(1,0)
        rod = Rod(start,end,youngs_modulus=E, cross_area_initial=A)
        
        matrix = np.zeros((4,4))
        matrix[0,0] = E*A/1
        matrix[2,0] = -E*A/1
        matrix[0,2] = -E*A/1
        matrix[2,2] = E*A/1


        nptest.assert_array_equal(matrix, element_matrix(rod)[0])

    def test_element_matrix_rod_rotation(self):
        E = 50e7
        A = .1
        angle = 1.12 # rad
        start = Node(0,0)
        end = Node(1*np.cos(angle),np.sin(angle))
        rod = Rod(start,end,youngs_modulus=E, cross_area_initial=A)
        L = np.sqrt(1)

        matrix = np.zeros((4,4))
        matrix[0,0] = E*A/L
        matrix[2,0] = -E*A/L
        matrix[0,2] = -E*A/L
        matrix[2,2] = E*A/L

        rot_mat = np.zeros((4,4))
        rot_mat[0,0] = np.cos(angle)
        rot_mat[1,0] = -np.sin(angle)
        rot_mat[0,1] = np.sin(angle)
        rot_mat[1,1] = np.cos(angle)
        rot_mat[2,2] = np.cos(angle)
        rot_mat[2,3] = np.sin(angle)
        rot_mat[3,2] = -np.sin(angle)
        rot_mat[3,3] = np.cos(angle)

        nptest.assert_array_equal(rot_mat, element_matrix(rod)[1])



    def test_solver(self):
        mesh5 = Mesh(0,0,0,0,initial_rod_area=.01, youngs_modulus=70e9)
        node_coords = [(0,0),(1,0)]
        node_connections = [(0,1)]
        
        mesh5.manual_mesh(node_coords, node_connections, seed = 1)

        support = Support((0,0),'p',orientation=0)
        mesh5.add_support(support)
        support = Support((1,0),'r',orientation=90)
        mesh5.add_support(support)

        force = Force(10000,0,(1,0))
        mesh5.add_external_force(force)

        displacements, forces = mesh5.solve(solver='rod', display_out=False)
        displacements_hand = np.zeros((4,1))
        displacements_hand[2,0] = 1.429e-5
        forces_hand = np.zeros((4,1))
        forces_hand[0,0] = -10000
        forces_hand[2,0] = 10000
        
        nptest.assert_array_almost_equal(displacements,displacements_hand)
        nptest.assert_array_almost_equal(forces,forces_hand)

        mesh6 = Mesh(0,0,0,0,initial_rod_area=.01, youngs_modulus=70e9)
        node_coords = [(0,0),(1,0),(0,1)]
        node_connections = [(0,1),(1,2),(0,2)]
        
        mesh6.manual_mesh(node_coords, node_connections, seed = 1)

        support = Support((0,0),'p',orientation=0)
        mesh6.add_support(support)
        support = Support((0,1),'r',orientation=0)
        mesh6.add_support(support)

        force = Force(10000,90,(1,0))
        mesh6.add_external_force(force)

        displacements, forces = mesh6.solve(solver='rod', display_out=False)
        displacements_hand = np.zeros((6,1))
        displacements_hand[2,0] = 1.429e-5
        displacements_hand[3,0] = 6.898e-5
        displacements_hand[5,0] = 1.429e-5
        forces_hand = np.zeros((6,1))
        forces_hand[0,0] = -10000
        forces_hand[1,0] = -10000
        forces_hand[3,0] = 10000
        forces_hand[4,0] = 10000
        
        nptest.assert_array_almost_equal(displacements,displacements_hand)
        nptest.assert_array_almost_equal(forces,forces_hand)
        

    def test_element_matrix_beam(self):
        E = 50e7
        A = .1
        I = .003
        L = 2
        start = Node(0,0)
        end = Node(L,0)
        rod = Rod(start,end,youngs_modulus=E, cross_area_initial=A, I = I)
        
        matrix = np.zeros((6,6))
        matrix[0,0] = E*A/L
        matrix[3,0] = -E*A/L
        matrix[0,3] = -E*A/L
        matrix[3,3] = E*A/L
        matrix[1,1] = 12*E*I/L**3
        matrix[1,2] = 6*E*I/L**2
        matrix[1,4] = -12*E*I/L**3
        matrix[1,5] = 6*E*I/L**2
        matrix[2,1] = 6*E*I/L**2
        matrix[2,2] = 4*E*I/L**1
        matrix[2,4] = -6*E*I/L**2
        matrix[2,5] = 2*E*I/L**1
        matrix[4,1] = -12*E*I/L**3
        matrix[4,2] = -6*E*I/L**2
        matrix[4,4] = 12*E*I/L**3
        matrix[4,5] = -6*E*I/L**2
        matrix[5,1] = 6*E*I/L**2
        matrix[5,2] = 2*E*I/L**1
        matrix[5,4] = -6*E*I/L**2
        matrix[5,5] = 4*E*I/L**1

        nptest.assert_array_equal(matrix, element_matrix_beam(rod)[0])

    def test_element_matrix_beam_rotated(self):
        E = 50e7
        A = .1
        I = .003
        L = 2
        start = Node(0,0)
        end = Node(0,L)
        rod = Rod(start,end,youngs_modulus=E, cross_area_initial=A, I = I)
        
        matrix = np.zeros((6,6))
        matrix[0,0] = E*A/L
        matrix[3,0] = -E*A/L
        matrix[0,3] = -E*A/L
        matrix[3,3] = E*A/L
        matrix[1,1] = 12*E*I/L**3
        matrix[1,2] = 6*E*I/L**2
        matrix[1,4] = -12*E*I/L**3
        matrix[1,5] = 6*E*I/L**2
        matrix[2,1] = 6*E*I/L**2
        matrix[2,2] = 4*E*I/L**1
        matrix[2,4] = -6*E*I/L**2
        matrix[2,5] = 2*E*I/L**1
        matrix[4,1] = -12*E*I/L**3
        matrix[4,2] = -6*E*I/L**2
        matrix[4,4] = 12*E*I/L**3
        matrix[4,5] = -6*E*I/L**2
        matrix[5,1] = 6*E*I/L**2
        matrix[5,2] = 2*E*I/L**1
        matrix[5,4] = -6*E*I/L**2
        matrix[5,5] = 4*E*I/L**1
        lam = 0
        mu  = 1
        T = np.eye(6)
        T[0,0]= lam
        T[1,1]= lam
        T[3,3]= lam
        T[4,4]= lam
        T[1,0]= -mu
        T[0,1]= mu
        T[4,3]= -mu
        T[3,4]= mu

        rotated_matrix = np.matmul(np.transpose(T),np.matmul(matrix,T))

        nptest.assert_array_equal(T, element_matrix_beam(rod)[1])

        nptest.assert_array_equal(rotated_matrix, element_matrix_beam(rod)[0])

if __name__ == '__main__':
    unittest.main()