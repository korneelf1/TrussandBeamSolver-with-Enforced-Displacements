# Start with rods, work towards beams

# Modules
import numpy as np
# import scipy as sp 
import matplotlib.pyplot as plt
import matplotlib.patches as pt
import math as m
import numpy.linalg as LA
import warnings
# initialize matplotlib
plt.figure()

# now tried to add rod generation in nodes part and gives singular matrix, check backupshit to see
class Node:
    '''Node object'''
    def __init__(self, x_coord, y_coord, forces_lst = None, support_lst =  None, moments_lst = None, adj_nodes = None, displacement = (0,0), inner = False):
        '''Creates an instance of a node and specifies connecting nodes.'''
        self.coord = (x_coord , y_coord)
        # self.adj_dict = {}
        self.connected_rods = []

        self.forces   = [] if forces_lst is None else forces_lst
        self.supports = [] if support_lst is None else support_lst
        self.moments  = [] if moments_lst is None else moments_lst
        self.internal_loads = [0,0,0]
        self.adj_nodes = {} if adj_nodes is None else adj_nodes
        self.displacement = displacement
        self.inner = inner # used for plotting these nodes smaller if necessary
        # for innnernodes
        self.main_node_start = None
        self.main_node_end = None

    def __str__(self) -> str:
        return str(self.coord)

    def neighbors(self):
        return [self.l_t,self.l, self.l_b, self.t, self.b, self.r_t, self.r, self.r_b]

    def internal_loads_calc(self):
        if self.main_node_start is not None:
            start = self.main_node_start
            
            self.internal_loads[2] = start.internal_loads[2] + (start.coord[0]-self.coord[0])*start.internal_loads[1] + (start.coord[1]-self.coord[1])*start.internal_loads[0] 

class Rod:
    def __init__(self, start_node, end_node, youngs_modulus = 70e9, cross_area_initial = .031416, min_area = 0, stress = 0, max_stress = 1e6, I = 0.0000785398, max_y = .1):
        '''Creates an instance of a rod between two nodes and specifies its' minimum area. Area in m^2.'''
        self.start = start_node
        self.end   = end_node
        self.__A_i   = cross_area_initial
        self.min_area = min_area
        self.E = youngs_modulus
        self.stress = stress
        self.stress_allowable = max_stress
        self.strain = 0
        self.I = I
        self.innerRods = []
        self.max_y = max_y
        
    def __str__(self) -> str:
        start = 'Start Node: ' + str(self.start) + '\n'
        end = 'End Node: ' + str(self.end) + '\n'
        length = 'Length: ' + str(self.length) + '\n'
        cross = 'Cross Section: ' + str(self.__A_i) 

        return str(start+end+length+cross)
    @property
    def A_i(self):
        '''Returns cross-sectional area of the rod.'''
        return self.__A_i
    
    @A_i.setter
    def A_i(self, value):
        '''Sets cross-sectional area of the rod if area greater than specified minimum'''
        self.__A_i = max(value, self.min_area)

    @property
    def length(self):
        return np.sqrt((self.start.coord[0] - self.end.coord[0])**2 + (self.start.coord[1] - self.end.coord[1])**2)

                
# start with rectangle
class Mesh:
    def __init__(self, lower_left, upper_left, lower_right, upper_right, max_displ = .01, initial_rod_area = .031416, youngs_modulus = 70e9, inertia =0.0000785398, max_y = .1):
        '''Enter corner nodes as tuples with x and y coordinate'''
        # self.ll_x, self.ll_y = lower_left[0], lower_left[1]
        # self.ul_x, self.ul_y = upper_left[0], upper_left[1]
        # self.lr_x, self.lr_y = lower_right[0], lower_right[1]
        # self.ur_x, self.ur_y = upper_right[0], upper_right[1]
        self.node_arr        = None
        self.displaced_nodes = []
        self.rod_lst         = []
        self.support_lst     = []
        self.forces_lst      = []
        self.moments_lst     = []
        self.nodes_dict = {}
        self.adj_dict = {}
        self.left_lower_node = None
        self.cross_area = initial_rod_area
        self.E = youngs_modulus
        self.__max_area = 0
        self.max_displ = max_displ
        self.inertia = inertia
        self.max_d = max_y

    def __repr__(self):
        return self.node_arr

    def manual_mesh(self, node_coords, node_connections, seed = 10):
        '''Enter node_coords as ordered list with (x,y) coordinates.
        Enter node_connections as list following order of node_coords with indices of connecting nodes.'''
        self.node_arr = []
        self.node_connections = node_connections

        self.seed = seed
        
        for i, node in enumerate(node_coords):
            new_node = Node(node[0], node[1])
            self.node_arr.append(new_node)

        for rod in node_connections:
            x_coord = (self.node_arr[rod[0]].coord[0],self.node_arr[rod[1]].coord[0])
            y_coord = (self.node_arr[rod[0]].coord[1],self.node_arr[rod[1]].coord[1])
            if seed == 1:
                new_rod = Rod(self.node_arr[rod[0]], self.node_arr[rod[1]],youngs_modulus= self.E,cross_area_initial=self.cross_area, I = self.inertia, max_y=self.max_d)
                self.node_arr[rod[0]].adj_nodes[str(self.node_arr[rod[1]].coord)] = new_rod
                self.node_arr[rod[1]].adj_nodes[str(self.node_arr[rod[0]].coord)] = new_rod
                

                self.rod_lst.append(new_rod)

                self.node_arr = np.array(self.node_arr)
            else:
                # beam
                x_coord_inter = np.linspace(x_coord[0],x_coord[1],seed+1)
                y_coord_inter = np.linspace(y_coord[0],y_coord[1],seed+1)
                
                # make all nodes, innernodes
                for i in range(len(x_coord_inter)):
                    if self.search_node((x_coord_inter[i], y_coord_inter[i])) is None:
                        new_node = Node(x_coord_inter[i], y_coord_inter[i], inner=True)
                        self.node_arr.append(new_node)
                        new_node.main_node_start = self.search_node((x_coord[0],y_coord[0]))
                        new_node.main_node_end = self.search_node((x_coord[1],y_coord[1]))


                # connect nodes
                for i in range(len(x_coord_inter)-1):
                    start_node = self.search_node((x_coord_inter[i], y_coord_inter[i]))
                    end_node = self.search_node((x_coord_inter[i+1], y_coord_inter[i+1]))                    

                    new_rod = Rod(start_node,end_node,youngs_modulus= self.E,cross_area_initial=self.cross_area, I = self.inertia, max_y=self.max_d)
                    self.rod_lst.append(new_rod)
                    start_node.adj_nodes[str((x_coord_inter[i+1], y_coord_inter[i+1]))] = new_rod
                    end_node.adj_nodes[str((x_coord_inter[i], y_coord_inter[i]))] = new_rod

        self.node_arr = np.array(self.node_arr)

    def search_node(self, coordinates, debug = False):
            '''Searches for node at given coordinates (tuple(x,y)). Returns None if not found.'''
            found_node = None
            
            if isinstance(self.node_arr, list):
                for node in self.node_arr:
                    if node.coord == coordinates:
                        found_node = node
                        if debug == True:
                            print(coordinates, node.coord)
                        return found_node
            else:
                for node in self.node_arr.flatten():
                    if node.coord == coordinates:
                        found_node = node
                        if debug == True:
                            print(coordinates, node.coord)
                        return found_node
                
            return found_node

            

    def generate_nodes(self, node_dist_hor, node_dist_vert):

        dist_h = node_dist_hor
        dist_v = node_dist_vert

        self.dx = dist_h
        self.dy = dist_v

        nodes_x = int((self.lr_x - self.ll_x) // node_dist_hor + 1)+1
        nodes_y = int((self.ul_y - self.ll_y) // node_dist_vert + 1)+1
        self.node_arr = np.full((nodes_x,nodes_y), None)
        self.coord_arr = np.zeros((nodes_x, nodes_y))
        # create lower_left node, from here every node can be accessed
       
        # adj_dict = {}
        left_n = None
        left_b_n = None
        bottom_n = None
        # slow algorithm
        # idea: set up random nodes in area connect them based on distance
        # could be totally random point locations or we could do rectengular pattern
        for j in range(nodes_y):
            # first traverse in x direction
            y_coord = self.ll_y + j * dist_v
            for i in range(nodes_x):
                x_coord = self.ll_x + i * dist_h
                
                
                curr = Node(x_coord, y_coord)
                if j == 0:
                    bottom_n = None
                else:
                    bottom_n = self.node_arr[j-1][i]
                    #stress initially 0 on default
                    rod = Rod(curr, bottom_n,youngs_modulus= self.E,cross_area_initial=self.cross_area)
                    bottom_n.adj_nodes[str(curr)] = rod
                    curr.adj_nodes[str(bottom_n)] = rod
                    self.rod_lst.append(rod)
                    # adj_dict[str(bottom_n)] = [bottom_n, {'area': self.cross_area,'E': self.E, 'stress': 0}]
                if i == 0:
                    left_n = None
                else:
                    left_n = self.node_arr[j][i-1]
                    rod = Rod(curr, left_n,youngs_modulus= self.E,cross_area_initial=self.cross_area)
                    left_n.adj_nodes[str(curr)] = rod
                    curr.adj_nodes[str(left_n)] = rod
                    self.rod_lst.append(rod)
                    # adj_dict[str(left_n)] = [left_n, {'area': self.cross_area,'E': self.E, 'stress': 0}]
                if i == 0 or j == 0:
                    left_b_n = None
                else: 
                    left_b_n = self.node_arr[j-1][i-1]
                    rod = Rod(curr, left_b_n,youngs_modulus= self.E,cross_area_initial=self.cross_area)
                    left_b_n.adj_nodes[str(curr)] = rod
                    curr.adj_nodes[str(left_b_n)] = rod
                    self.rod_lst.append(rod)
                    # adj_dict[str(left_b_n)] = [left_b_n, {'area': self.cross_area,'E': self.E, 'stress': 0}]
                

                curr.left = left_n
                curr.left_bottom = bottom_n

                self.node_arr[j][i] = curr

                # to make graphs 
                self.nodes_dict[str(x_coord) + ',' + str(y_coord)] = curr
                # create adjacent edges

        self.left_lower_node = self.node_arr[0][0]
        self.node_arr = self.node_arr.flatten()

    def generate_rods(self,min_area = 0):
        node = self.node_arr[0][0]
        node_row = 0
        # until we are in rightmost upper node
        end = False
        while not end:
            coord_curr = node.coord
            if node.r is not None:
                rod = Rod(self.search_node((coord_curr)), self.search_node((node.r.coord)))

                self.rod_lst.append(rod)
                node_next = node.r
                if node.t is not None:
                    rod = Rod(self.search_node((node.r.coord)), self.search_node((node.t.coord)))
                    self.rod_lst.append(rod)                   
            else:
                node_row += 1
                if node_row < len(self.node_arr):
                    node_next = self.node_arr[node_row][0]
                else: 
                    end = True
            if node.t is not None:
                rod = Rod(self.search_node((coord_curr)), self.search_node((node.t.coord)))
                self.rod_lst.append(rod)

            if node.r_t is not None:
                rod = Rod(self.search_node((coord_curr)), self.search_node((node.r_t.coord)))
                self.rod_lst.append(rod)

            node = node_next
            

    def add_external_force(self, force):
        '''Adds external force, enter Force object.'''
        
        # drawing 
        force_coord = force.node
        length      = force.magn / 10000
        dx = m.cos(force.angle)* length
        dy = m.sin(force.angle) * length
        # plt.arrow(force_coord[0], force_coord[1], dx, dy, color ='r', width = .015)
        plt.text(force_coord[0],force_coord[1],force.magn)
        # search and append force to node
        node_force = self.search_node(force.node)
        node_force.forces.append(force)

        # append force to mesh force list
        self.forces_lst.append(force)


    def add_external_moment(self, moment):
        '''Adds external moment, enter Moment object.'''
        # drawing 
        mom_coord = moment.node
        
        # color 
        c = 'r'
        if moment.magn < 0:
            c = 'k'
        circle = plt.Circle(mom_coord, radius = 1, color=c,fill = False)
        # retrieve axes to plot circle
        ax = plt.gca()
        ax.add_patch(circle)
        # search and append force to node
        node_moment = self.search_node(mom_coord)
        node_moment.moments.append(moment)

        # append force to mesh force list
        self.moments_lst.append(moment)
        

    def add_support(self, support):
        '''Adds support, enter Support object.'''
        if support.valid:
            # drawing
            con_coord = support.node
            plt.scatter(con_coord[0], con_coord[1], c='b')
            support_node = self.search_node(support.node)
            if support_node:
                support_node.supports.append(support)
                self.support_lst.append(support)

            else:
                warnings.warn('The support is not connected to the structure.')
                
        else:
            warnings.warn('Support is not valid.')

    def solve(self, solver= 'beam', display_out = False, target_displ = False):
        '''Solves Rods given forces, moments and supports.
        Returns None if not solveable'''
        # delete node if rods are all gone

        # def element_matrix(rod):
        #     '''Returns 4x4 matrix in global coordinates.'''
        #     E = self.E
        #     rot_mat = np.zeros((2,2))
        #     start_node = rod.start.coord
        #     end_node   = rod.end.coord

        #     U1 = np.array([end_node[0] - start_node[0],end_node[1] - start_node[1]])
        #     U1_norm = np.sqrt(np.sum(U1**2)) # this is length of the element
        #     U1 = U1/U1_norm

        #     U2 = np.array([-U1[1], U1[0]])

        #     rot_mat[0], rot_mat[1] = U1, U2
        #     rot_mat = rot_mat.transpose() # why here suddenly rotated

        #     local_K = np.zeros((4,4))
        #     local_K[0][0] = E*rod.A_i/U1_norm
        #     local_K[0][2] = -E*rod.A_i/U1_norm
        #     local_K[2][0] = -E*rod.A_i/U1_norm
        #     local_K[2][2] = E*rod.A_i/U1_norm

        #     T1 = np.zeros((4,4))
        #     T1[0:2,0:2] = rot_mat
        #     T1[2:,2:] = rot_mat
        #     T2 = np.zeros((4,4))
        #     T2[0:2,0:2] = rot_mat.transpose()
        #     T2[2:,2:] = rot_mat.transpose()

        #     return np.matmul(T1,np.matmul(local_K,T2))

        def element_matrix(rod):
            '''Returns 4x4 matrix in global coordinates.'''
            E = rod.E
            rot_mat = np.zeros((2,2))
            start_node = rod.start.coord
            end_node   = rod.end.coord

            U1 = np.array([end_node[0] - start_node[0],end_node[1] - start_node[1]])
            U1_norm = np.sqrt(np.sum(U1**2)) # this is length of the element
            U1 = U1/U1_norm

            U2 = np.array([-U1[1], U1[0]])

            rot_mat[0], rot_mat[1] = U1, U2
            # rot_mat = rot_mat.transpose()

            local_K = np.zeros((4,4))
            local_K[0][0] = E*rod.A_i/U1_norm
            local_K[0][2] = -E*rod.A_i/U1_norm
            local_K[2][0] = -E*rod.A_i/U1_norm
            local_K[2][2] = E*rod.A_i/U1_norm

            T1 = np.zeros((4,4))
            T1[0:2,0:2] = rot_mat
            T1[2:,2:] = rot_mat
    

            return np.matmul(np.transpose(T1),np.matmul(local_K,T1))

        def element_matrix_beam(beam):
            E = self.E
            A = beam.A_i
            I = beam.I
            rot_mat = np.zeros((3,3))
            start_node = beam.start.coord
            end_node   = beam.end.coord

            U1 = np.array([end_node[0] - start_node[0],end_node[1] - start_node[1]])
            U1_norm = np.sqrt(np.sum(U1**2)) # this is length of the element
            U1 = U1/U1_norm

            U2 = np.array([-U1[1], U1[0]])

            rot_mat[0][0:2], rot_mat[1][0:2] = U1, U2
            # rot_mat = rot_mat.transpose()
            rot_mat[2,2] = 1

            T1 = np.zeros((6,6))
            T1[:3,:3] = rot_mat
            T1[3:,3:] = rot_mat


            local_K = np.zeros((6,6))
            local_K[0,0], local_K[3,3] = E*A/U1_norm, E*A/U1_norm
            local_K[0,3], local_K[3,0] = -E*A/U1_norm, -E*A/U1_norm
            local_K[1,1], local_K[4,4] = 12*E*I/U1_norm**3,12*E*I/U1_norm**3
            local_K[1,4], local_K[4,1] = -12*E*I/U1_norm**3,-12*E*I/U1_norm**3
            local_K[1,2], local_K[2,1], local_K[1,5], local_K[5,1] = 6*E*I/U1_norm**2, 6*E*I/U1_norm**2, 6*E*I/U1_norm**2, 6*E*I/U1_norm**2
            local_K[2,4], local_K[4,2], local_K[4,5], local_K[5,4] = -6*E*I/U1_norm**2, -6*E*I/U1_norm**2, -6*E*I/U1_norm**2, -6*E*I/U1_norm**2
            local_K[2,2], local_K[5,5] = 4*E*I/U1_norm, 4*E*I/U1_norm    
            local_K[2,5], local_K[5,2] = 2*E*I/U1_norm, 2*E*I/U1_norm

            global_K = np.matmul(np.transpose(T1),np.matmul(local_K,T1))

            return global_K


        def internal_stress(node_array, displacements, forces, solver = 'rod'):
            '''Innefficient part: all rods cylced over twice'''
            rods_passed = []
            for i, node in enumerate(node_array):
                for rod in node.adj_nodes.values():
                    if rod not in rods_passed:
                        if solver =='rod':
                            start_disp = (displacements[2*i],displacements[2*i+1])
                            j = np.where(node_array == rod.end)[0]
                            end_disp = (displacements[2*j],displacements[2*j+1])

                            tot_elongation = np.sqrt((end_disp[0] - start_disp[0])**2 + (end_disp[1] - start_disp[1])**2)
                            original_length = rod.length
                            strain = tot_elongation/original_length
                            # strain = (tot_elongation-original_length)/original_length

                            rod.strain = strain
                            sigma = rod.E*strain
                            rod.stress = sigma[0][0]
                            rods_passed.append(rod)
                        else:
                            start_disp = (displacements[3*i],displacements[3*i+1])
                            j = np.where(node_array == rod.end)[0]
                            end_disp = (displacements[3*j],displacements[3*j+1])

                            tot_elongation = np.sqrt((end_disp[0] - start_disp[0])**2 + (end_disp[1] - start_disp[1])**2)
                            original_length = rod.length
                            strain = (tot_elongation)/original_length
                            rod.strain = strain
                            sigma = rod.E*strain
                            # rod.stress = sigma[0][0]
                            rods_passed.append(rod)
                            # bending stress
                            node.internal_loads_calc()
                            mom = rod.start.internal_loads[2]
                            bending_stress = mom/rod.I*rod.max_y

                            if abs(bending_stress + rod.stress)>(-bending_stress + rod.stress):
                                rod.stress = bending_stress + rod.stress
                            else:
                                rod.stress = -bending_stress + rod.stress
    	                    

        N_rods = len(self.rod_lst)
        N_nodes = len(self.node_arr.flatten())
        # flatten node_arr
        node_vector = self.node_arr.reshape(N_nodes)
            

        if solver == 'rod':
            # every node has 2 DoF (in case of two dims)
            K = np.zeros((2*N_nodes,2*N_nodes)) # stiffness matrix
            u = np.ones((2*N_nodes,1)) # displacement vector, initialize with ones, where constrained assign 0, now we find constraint DoF easily
            f = np.zeros((2*N_nodes,1)) # force vector

            
            # loop over rods
            for i_rod in range(N_rods):
                curr_rod = self.rod_lst[i_rod]

                K_rod = element_matrix(curr_rod)
                
                # assembling in the global K matrix
                start_node = curr_rod.start
                end_node  = curr_rod.end

                ind_start_node = np.where(node_vector == start_node)[0][0]
                ind_end_node = np.where(node_vector == end_node)[0][0]

                for i, index_x in enumerate([2*ind_start_node,2*ind_start_node+1, 2*ind_end_node,2*ind_end_node+1]):
                    for j, index_y in enumerate([2*ind_start_node,2*ind_start_node+1, 2*ind_end_node,2*ind_end_node+1]):       
                        K[index_y,index_x] += K_rod[j,i]

            # apply boundary conditions
            for support in self.support_lst:
                support_node = self.search_node(support.node)
                ind_support_node = np.where(node_vector == support_node)[0][0]
                if support.orientation%(np.pi) == 0: # along x - axis
                    u[2*ind_support_node] = 0
                    if support.par:
                        u[2*ind_support_node + 1] = 0
                elif support.orientation%(np.pi/2) == 0: # along y - axis
                    u[2*ind_support_node + 1] = 0
                    if support.par:
                        u[2*ind_support_node] = 0

            if target_displ is not False:
                for node_coord, target in target_displ.items():
                    displ_node = self.search_node(eval(node_coord))

                    ind_displ_node = np.where(node_vector == displ_node)[0]
                    
                    u[2*ind_displ_node] = target[0]
                    u[2*ind_displ_node + 1] = target[1]
                            
            print('line 574')
            if len(u[u == 0]) <3:
                print('Structure is underconstrained.')
            # apply external forces
            for force in self.forces_lst:
                force_node = self.search_node(force.node)
                ind_force_node = np.where(node_vector == force_node)[0][0]
                fx = force.magn * np.cos(force.angle)
                fy = force.magn * np.sin(force.angle)

                f[2*ind_force_node] = fx
                f[2*ind_force_node + 1] = fy


            K_alt = np.copy(K)
            # rearranging matrices making reduced matrix
            ind_bcs = np.where(u!=1)[0]
            ind_unkn = np.where(u==1)[0]
            K_red = K_alt[ind_unkn,:]
            u_red = u[ind_bcs]
            K_red = K_red[:,ind_bcs]
            correction = np.matmul(K_red,u_red)
            for i, index in enumerate(ind_unkn):
                # print(index,i)
                f[index] -= correction[i][0]
            # f_red = f[ind_help]




            for i in ind_bcs:
                new_row = np.zeros(len(K[0]))
                new_row[i] = 1
                K_alt[i] = new_row
                f[i] = 0
            

            displ = LA.solve(K_alt,f)

            if target_displ is not False:
                for node_coord, target in target_displ.items():
                    displ_node = self.search_node(eval(node_coord))

                    ind_displ_node = np.where(node_vector == displ_node)[0]
                    
                    displ[2*ind_displ_node] = target[0]
                    displ[2*ind_displ_node + 1] = target[1]
            forces = np.matmul(K,displ)


            

        elif solver == 'beam':
            # every node has 3 DoF (in case of three dims)
            # nodes including innernodes
            els = 3*(N_nodes)
            K = np.zeros((els,els)) # stiffness matrix
            u = np.ones((els,1)) # displacement vector, initialize with ones, where constrained assign 0, now we find constraint DoF easily
            f = np.zeros((els,1)) # force and moment vector

            # loop over rods
            for i_rod in range(N_rods):
                curr_rod = self.rod_lst[i_rod]

                
                K_rod = element_matrix_beam(curr_rod)
                # assembling in the global K matrix
                start_node = curr_rod.start 
                end_node  = curr_rod.end 

                ind_start_node = np.where(node_vector == start_node)[0][0]
                ind_end_node = np.where(node_vector == end_node)[0][0]

                for i, index_x in enumerate([3*ind_start_node,3*ind_start_node+1,3*ind_start_node+2, 3*ind_end_node,3*ind_end_node+1,3*ind_end_node+2]):
                    for j, index_y in enumerate([3*ind_start_node,3*ind_start_node+1,3*ind_start_node+2, 3*ind_end_node,3*ind_end_node+1,3*ind_end_node+2]):       
                        K[index_y,index_x] += K_rod[j,i]


            # apply boundary conditions
            for support in self.support_lst:
                support_node = self.search_node(support.node)
                ind_support_node = np.where(node_vector == support_node)[0][0]
                if support.orientation%(np.pi) == 0: # along x - axis
                    u[3*ind_support_node] = 0
                    if support.par:
                        u[3*ind_support_node + 1] = 0
                elif support.orientation%(np.pi/2) == 0: # along y - axis
                    u[3*ind_support_node + 1] = 0
                    if support.par:
                        u[3*ind_support_node] = 0
                if support.type == 'f':
                    u[3*ind_support_node + 2] = 0

            if len(u[u == 0]) <3:
                print('Structure is underconstrained.')
            # apply external forces
            for force in self.forces_lst:
                force_node = self.search_node(force.node)
                ind_force_node = np.where(node_vector == force_node)[0][0]
                fx = force.magn * np.cos(force.angle)
                fy = force.magn * np.sin(force.angle)

                f[3*ind_force_node] = fx
                f[3*ind_force_node + 1] = fy

            # apply external moments
            for moment in self.moments_lst:
                # print(moment.node)
                moment_node = self.search_node(moment.node)
                ind_moment_node = np.where(node_vector == moment_node)[0][0]
                # print(ind_moment_node)
                f[3*ind_moment_node + 2] = moment.magn
            
            K_alt = np.copy(K)

            # # rearranging matrices
            # ind_bcs = np.where(u==0)[0]
            if target_displ is not False:
                for node_coord, target in target_displ.items():
                    displ_node = self.search_node(eval(node_coord))

                    ind_displ_node = np.where(node_vector == displ_node)[0]
                    
                    u[3*ind_displ_node] = target[0]
                    u[3*ind_displ_node + 1] = target[1]

            # rearranging matrices making reduced matrix
            
            ind_bcs = np.where(u!=1)[0]
            ind_unkn = np.where(u==1)[0]
            K_red = K_alt[ind_unkn,:]
            u_red = u[ind_bcs]
            K_red = K_red[:,ind_bcs]
            correction = np.matmul(K_red,u_red)
            for i, index in enumerate(ind_unkn):
                f[index] -= correction[i][0]

            for i in ind_bcs:
                new_row = np.zeros(len(K[0]))
                new_row[i] = 1
                K_alt[i] = new_row
                f[i] = 0

            displ = LA.solve(K_alt,f)
            if target_displ is not False:
                for node_coord, target in target_displ.items():
                    displ_node = self.search_node(eval(node_coord))

                    ind_displ_node = np.where(node_vector == displ_node)[0]
                    
                    displ[3*ind_displ_node] = target[0]
                    displ[3*ind_displ_node + 1] = target[1]
            
            forces = np.matmul(K,displ)
            
        for i, node in enumerate(node_vector):
            index_corr	= 2
            internal_loads = [forces[index_corr*i], forces[index_corr*i+1],0]
            if solver == 'beam':
                index_corr = 3
                internal_loads = [forces[index_corr*i], forces[index_corr*i+1],forces[index_corr*i+2]]

            displacement = (displ[index_corr*i],displ[index_corr*i+1])
            node.displacement = displacement
            
            node.internal_loads = internal_loads
            node.internal_loads_calc()

        internal_stress(node_vector, displ, forces, solver=solver)

        if display_out:
            print('Displacements: ')
            # print(displ)
            for node in self.node_arr:
                if node.inner is False:
                    print(node.displacement)
            print('\nForces: ')
            print(forces)
            print('\nRods strains: ')
            for rod in self.rod_lst:
                print(rod.strain)
            print('\nRods stresses: ')
            for rod in self.rod_lst:
                print(rod.stress) 
            print('\nRods areas: ')
            for rod in self.rod_lst:
                print(rod.A_i)
        stresses = []
        for rod in self.rod_lst:
            # print(rod)
            max_stress = rod.stress_allowable
            stresses.append(rod.stress)
        if np.max(np.abs(stresses)) > max_stress:
            print('The maximum stress exceeds the allowable stress.')
        return displ, forces


    def optimize(self, optimizer = 'rod', max_iter = 15, area_treshold = 1e-8, stress_treshold = 1e-5, target_displ = False):
        # current fix
        max_displ = self.max_displ
        '''Area_treshold is minimum area before deleted'''

        if target_displ is False:
            for i in range(max_iter):
                self.solve(solver = optimizer,target_displ=target_displ)

                print('max area is in!')
                max_area = 1
                self.__max_area = 0 # used to scale values for transparancy in plot            
                for rod in self.rod_lst:
                    if rod.A_i > area_treshold:
                        new_area = rod.stress/rod.stress_allowable*rod.A_i
                        if new_area > max_area:
                            new_area = 1
                        if new_area > area_treshold:
                            rod.A_i = new_area
                            self.__max_area = max(new_area, self.__max_area)
                        else:
                            self.delete_rod(rod)
                    else:
                        self.delete_rod(rod)

                    
        else:
            for i in range(max_iter):
                self.solve(solver = optimizer,target_displ=target_displ)
                for rod in self.rod_lst:
                    if rod.A_i > area_treshold:
                        new_area = rod.stress/rod.stress_allowable*rod.A_i
                        
                        rod.A_i = new_area
                      
                    else:
                        self.delete_rod(rod)
            


    def delete_rod(self, rod):
        # problem: rods are generated and put twice in the rod_lst 
        self.rod_lst.remove(rod)
        print(rod.start.adj_nodes)
        del rod.start.adj_nodes[str(rod.end)]
        if not rod.start.adj_nodes:
            # dict is empty no rods connecting node
            self.node_arr = self.node_arr[self.node_arr!= rod.start]

    
        del rod.end.adj_nodes[str(rod.start)]
        if not rod.end.adj_nodes:
            # dict is empty no rods connecting node
            self.node_arr = self.node_arr[self.node_arr!= rod.end]

    def plot(self, displ_factor = 1, plot_inner = False, arrow_width = .015, arrow_length = 100000):
        '''Plots mesh, displacement factor scales the displacement for visual purposes'''
        # plotting rods initial
        for rod in self.rod_lst:
            transparancy = max(.1, rod.A_i/max(rod.A_i,self.__max_area))
            plt.plot((rod.start.coord[0],rod.end.coord[0]),(rod.start.coord[1],rod.end.coord[1]), 'k', alpha = transparancy)
            x_coord = (rod.start.coord[0] + rod.start.displacement[0]*displ_factor,rod.end.coord[0] + rod.end.displacement[0]*displ_factor)
            y_coord = (rod.start.coord[1] + rod.start.displacement[1]*displ_factor,rod.end.coord[1] + rod.end.displacement[1]*displ_factor)

            plt.plot(x_coord,y_coord, c='c', alpha = transparancy, label = rod.stress)

        for node in self.node_arr.flatten():
            if node.inner:
                if plot_inner:
                    # before loads are applied
                    plt.scatter(*node.coord, c='r', s = 25)
                    # location after loads are applied
                    x_coord = node.coord[0] + node.displacement[0]*displ_factor
                    y_coord = node.coord[1] + node.displacement[1]*displ_factor

                    plt.scatter(x_coord, y_coord, color = 'g', s = 10)
            
            else:
                
                # before loads are applied
                plt.scatter(*node.coord, c='y', s = 100)
                # location after loads are applied
                x_coord = node.coord[0] + node.displacement[0]*displ_factor
                y_coord = node.coord[1] + node.displacement[1]*displ_factor

                plt.scatter(x_coord, y_coord, color = 'g')
        for force in self.forces_lst:
            # drawing 
            force_coord = force.node
            length      = force.magn / arrow_length
            dx = m.cos(force.angle)* length
            dy = m.sin(force.angle) * length
            plt.arrow(force_coord[0], force_coord[1], dx, dy, color ='r', width = arrow_width)
            plt.text(force_coord[0],force_coord[1],force.magn)
        plt.gca().set_aspect('equal')
        
        plt.show()

    def plot_shear_moment(self):
        '''Plots mesh, displacement factor scales the displacement for visual purposes'''
        shear = []
        moments = []
        for node in self.node_arr:
            shear, moments = node.internal_loads[1][0], node.internal_loads[2][0]
            x,y = node.coord[0], node.coord[1]
            plt.scatter(x,y+shear/1,c='r')

            plt.scatter(x,y+moments/1,c='g')

        
        plt.show()

# for unit testing
def element_matrix(rod):
        '''Returns 4x4 matrix in global coordinates.'''
        E = rod.E
        rot_mat = np.zeros((2,2))
        start_node = rod.start.coord
        end_node   = rod.end.coord

        U1 = np.array([end_node[0] - start_node[0],end_node[1] - start_node[1]])
        U1_norm = np.sqrt(np.sum(U1**2)) # this is length of the element
        U1 = U1/U1_norm

        U2 = np.array([-U1[1], U1[0]])

        rot_mat[0], rot_mat[1] = U1, U2
        # rot_mat = rot_mat.transpose()

        local_K = np.zeros((4,4))
        local_K[0][0] = E*rod.A_i/U1_norm
        local_K[0][2] = -E*rod.A_i/U1_norm
        local_K[2][0] = -E*rod.A_i/U1_norm
        local_K[2][2] = E*rod.A_i/U1_norm

        T1 = np.zeros((4,4))
        T1[0:2,0:2] = rot_mat
        T1[2:,2:] = rot_mat
        print(T1)
        # T2 = np.zeros((4,4))
        # T2[0:2,0:2] = rot_mat.transpose()
        # T2[2:,2:] = rot_mat.transpose()

        return np.matmul(np.transpose(T1),np.matmul(local_K,T1)), T1

def element_matrix_beam(beam):
    E = beam.E
    A = beam.A_i
    I = beam.I
    rot_mat = np.zeros((3,3))
    start_node = beam.start.coord
    end_node   = beam.end.coord

    U1 = np.array([end_node[0] - start_node[0],end_node[1] - start_node[1]])
    U1_norm = np.sqrt(np.sum(U1**2)) # this is length of the element
    U1 = U1/U1_norm

    U2 = np.array([-U1[1], U1[0]])

    rot_mat[0][0:2], rot_mat[1][0:2] = U1, U2
    # rot_mat = rot_mat.transpose()
    rot_mat[2,2] = 1

    T1 = np.zeros((6,6))
    T1[:3,:3] = rot_mat
    T1[3:,3:] = rot_mat


    local_K = np.zeros((6,6))
    local_K[0,0], local_K[3,3] = E*A/U1_norm, E*A/U1_norm
    local_K[0,3], local_K[3,0] = -E*A/U1_norm, -E*A/U1_norm
    local_K[1,1], local_K[4,4] = 12*E*I/U1_norm**3,12*E*I/U1_norm**3
    local_K[1,4], local_K[4,1] = -12*E*I/U1_norm**3,-12*E*I/U1_norm**3
    local_K[1,2], local_K[2,1], local_K[1,5], local_K[5,1] = 6*E*I/U1_norm**2, 6*E*I/U1_norm**2, 6*E*I/U1_norm**2, 6*E*I/U1_norm**2
    local_K[2,4], local_K[4,2], local_K[4,5], local_K[5,4] = -6*E*I/U1_norm**2, -6*E*I/U1_norm**2, -6*E*I/U1_norm**2, -6*E*I/U1_norm**2
    local_K[2,2], local_K[5,5] = 4*E*I/U1_norm, 4*E*I/U1_norm    
    local_K[2,5], local_K[5,2] = 2*E*I/U1_norm, 2*E*I/U1_norm

    global_K = np.matmul(T1.transpose(),np.matmul(local_K,T1))

    return global_K, T1





if __name__ == '__main__':
    mesh = Mesh((0,0),(0,1),(1,0),(1,1))
    # mesh.generate_nodes(.25,.25)




