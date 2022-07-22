import MeshGenerator as mg
import math as m
import matplotlib.pyplot as plt
import warnings

class Force:
    def __init__(self, magnitude, orientation, application_node):
        '''Give magnitude in Newton, orientation in degrees wrt to horizontal x-axis. Application node as tuple of indices in node array'''
        self.magn = magnitude
        self.angle = m.radians(orientation)
        self.node  = application_node


    def draw(self, mesh):
        force_coord = mesh.node_arr[self.node[0]][self.node[1]].coord
        length      = self.magn / 100
        dx = m.cos(self.angle)* length
        dy = m.sin(self.angle) * length
        plt.arrow(force_coord[0], force_coord[1], dx, dy, color ='r')

    def __repr__(self) -> str:
        return str(f'magnitude = {self.magn}\norientation in rad = {self.angle}\nnode = {self.node}')

class Moment:
    '''Enter magnitude (pos is ccw) in Newtonmeter and application node.'''
    def __init__(self, magnitude, application_node):
        self.magn = magnitude
        self.node = application_node

    

class Support:
    '''Application node as tuple of indices.
    Types of supports:
    Roller ('r'): can support perpendicular force
    Fixed ('f'): can support perpendicular force, parallel force, moment
    Pinned ('p'): can support perpendicular force, parallel force
    Simple ('s'): can support perpendicular force, parallel force in case of friction 
    (Orientation in degrees, give angle of perpendicular force wrt horizontal x-axis)'''
    def __init__(self, application_node, type, orientation = 90):
        self.node = application_node
        self.type = type
        self.orientation = m.radians(orientation)

        self.valid = True

        self.perp = True
        if self.type in ('f', 'p', 's'):
            self.par = True
        else:
            self.par = False
        if self.type not in ('f', 'p', 's','r'):
            warnings.warn('Support type does not exist, try again.')
            self.valid=False

    def draw(self, mesh):
        con_coord = mesh.node_arr[self.node[0]][self.node[1]].coord
        plt.scatter(con_coord[0], con_coord[1], c='b')

if __name__ == "__main__":	
    
    '''
    '''
    def example_1():
        mesh = mg.Mesh((0,0),(0,10),(10,0),(10,10))
        # example 1
        # mesh.generate_nodes(5,5)
        # # mesh.generate_rods()

        force1 = Force(10000, 0, (0,10))
        # mesh.add_external_force(force1)
        # force2 = Force(10000, 270, (10,10))
        # mesh.add_external_force(force2)

        # support1 = Support((10.,5), 'f')
        # support2 = Support((10,0), 'f')
        # support3 = Support((0.,0), 'r', orientation=0)
        # mesh.add_support(support1)
        # mesh.add_support(support2)
        # mesh.add_support(support3)

    def example_rod_2():
        mesh = mg.Mesh((0,0),(0,10),(10,0),(10,10))
        node_coords = [(0,0),(1,0),(0,1),(1,1)]
        node_connections = [(0,1),(1,2),(0,2),(1,3),(2,3)]
        mesh.manual_mesh(node_coords, node_connections, seed = 1)

        force1 = Force(7.54618683e+03, 90, (1,1))
        # mesh.add_external_force(force1)


        support1 = Support((0,0),'p')
        support2 = Support((0,1),'r', orientation=0)
        mesh.add_support(support1)
        mesh.add_support(support2)

        mesh.solve(solver='rod', display_out=True,target_displ={'(1,1)':[0e-5,2e-5]})
        mesh.optimize(target_displ={'(1,1)':[0e-5,2e-5]})
        mesh.solve(solver='rod', display_out=True,target_displ={'(1,1)':[0e-5,2e-5]})

        # mesh.solve(solver='rod', display_out=True)

        mesh.plot(displ_factor=1000)
    
    def example_3():
        mesh = mg.Mesh((0,0),(0,10),(10,0),(10,10))
        node_coords = [(0,0),(1,0)]
        node_connections = [(0,1)]
        mesh.manual_mesh(node_coords, node_connections, seed=5)
        support1 = Support((0,0),'f')
        mesh.add_support(support1)

        force = Force(1.64933580e+04, 90, (1,0))
        mesh.add_external_force(force)
        # mesh.solve(solver = 'beam', display_out=True, target_displ={'(1,0)':[0,1e-3]})
        mesh.solve(solver = 'beam', display_out=True)
        # mesh.plot_shear_moment()
        mesh.plot(displ_factor=50)

    def validation_rod_1():
        '''
        Example from svv reader
        '''
        mesh = mg.Mesh((0,0),(0,.300),(.300,0),(.300,.300), youngs_modulus=210e9,initial_rod_area=400e-6)
        node_coords = [(0,0),(.100,0), (.200,0), (.100,.100)]
        node_connections = [(0,1), (0,3),(1,2),(1,3),(2,3)]
        mesh.manual_mesh(node_coords, node_connections, seed=1)

        support1 = Support((0,0),'p')
        mesh.add_support(support1)
        support2 = Support((.200,0),'p')
        mesh.add_support(support2)
        

        force = Force(1000, 270, (.100,0))
        mesh.add_external_force(force)

        mesh.solve(solver = 'rod', display_out=True)
        mesh.plot(arrow_width=.0015, displ_factor=1000)

    def example_beam_2():
        mesh = mg.Mesh((0,0),(100,0),(0,100),(100,100), youngs_modulus=200e9,initial_rod_area=400e-4)
        node_coords = [(0,0),(0, 10)]
        node_connections = [(0,1)]
        mesh.manual_mesh(node_coords, node_connections, seed=10)
        support1 = Support((0,0),'p')
        mesh.add_support(support1)
        support2 = Support((0,10),'r', 0)
        mesh.add_support(support2)

        force = Force(10000, 90, (10,0))

        moment = Moment(1000000,(0,0))
        mesh.add_external_moment(moment)
        mesh.solve(solver = 'beam', display_out=False)
        mesh.plot(displ_factor=10)

    def example_beam_3():
        mesh = mg.Mesh((0,0),(100,0),(0,100),(100,100), youngs_modulus=200e9,initial_rod_area=400e-4)
        node_coords = [(0,0),(0, 10),(10,0),(2,5),(10,10)]
        node_connections = [(0,1),(0,2),(1,2),(0,3),(2,3),(3,4),(4,1),(4,2)]
        mesh.manual_mesh(node_coords, node_connections, seed=20)
        support1 = Support((0,0),'p')
        mesh.add_support(support1)
        support2 = Support((0,10),'r', 0)
        mesh.add_support(support2)

        force = Force(10000, 270, (10,0))
        # mesh.add_external_force(force)

        moment = Moment(1000,(0,0))
        mesh.add_external_moment(moment)
        # mesh.plot()
        mesh.solve(solver = 'beam', display_out=False, target_displ={'(10,10)':[0e-5,0e-5]})
        mesh.optimize(optimizer = 'beam',max_iter = 10)
        mesh.solve(solver = 'beam', display_out=False, target_displ={'(10,10)':[0e-5,0e-5]})
        mesh.plot(displ_factor=10000, plot_inner=False)


    def validation_beam_1():
        mesh = mg.Mesh(0,0,0,0,inertia=1.333e-13,max_y=.0002, youngs_modulus=194.3e9, initial_rod_area=1e-5)
        node_coords = [(0,0),(.4,0)]
        node_connections = [(0,1)]
        mesh.manual_mesh(node_coords,node_connections, seed = 10)
        
        support1 = Support((0,0),'f')
        mesh.add_support(support1)

        force = Force(.196,orientation=270,application_node=(.4,0))
        mesh.add_external_force(force)

        mesh.solve(solver = 'beam', display_out=True)
        mesh.plot(displ_factor=1)

    def validation_beam_2():
        mesh = mg.Mesh(0,0,0,0,inertia=10000e-12,max_y=.0002, youngs_modulus=70e9, initial_rod_area=100e-6)
        node_coords = [(0,0),(.075,0),(.15,0)]
        node_connections = [(0,1),(1,2)]
        mesh.manual_mesh(node_coords,node_connections, seed =2)
        
        support1 = Support((0,0),'p')
        mesh.add_support(support1)
        support2 = Support((0.15,0),'p')
        mesh.add_support(support2)

        force = Force(1000,orientation=270,application_node=(.075,0))
        mesh.add_external_force(force)

        mesh.solve(solver = 'beam', display_out=True)
        mesh.plot(displ_factor=100, arrow_width=.0015)
        
    validation_beam_2()
    
    # mesh.optimize(max_iter=100)
    # mesh.solve(solver = 'rod', display_out=True)
    # mesh.plot()

       