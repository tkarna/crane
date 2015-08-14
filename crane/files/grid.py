"""
Grid objects

Method to read and write grid data.

Open boundary names are hard-coded here. If new open boundaries are added to a
grid they should be added here so the names will be added to the grid object.
genBCTides relies on having these names.

lopezj - 08/06/2012
"""

#-------------------------------------------------------------------------------
# Imports
#-------------------------------------------------------------------------------
import numpy as np

#-------------------------------------------------------------------------------
# Constants 
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Classes and functions
#-------------------------------------------------------------------------------
class Object(object):
    pass

def _readBndryNodes(f):
    """ 
        Reads the boundary section of hgrid.gr3 and returns nodes and 
        a information about the boundary type and name.

        obNodes - Object with a collection of boundaries
            ob_?_nodes - Array of nodes along boundary number ?
            ob_?_type - Boundary type {'ocean','river','land','island'}
            ob_?_name - Open boundary name 
    """
    #OB-Open boundary, code copy and pasted from there, but works generically
    _tmp = f.readline()
    # In case islands or land is missing, I guess you can have two land 
    # boundaries only, but I'll pretend that will not happen
    if _tmp == '':
        return
    nOB = int(_tmp.rstrip("\n").split()[0])
    nNodesTotalOB = int(f.readline().rstrip("\n").split()[0])
    obNodes = Object()
    obNodes.nB = nOB
    # Loop over all boundaries
    for ob in range(nOB):
        str = f.readline().rstrip("\n")
        nNodesThisOB = int(str.split()[0])
        # Identify as ocean, river, land, or island
        str = ' '.join(str.split()[2:])      
        if str.find('ocean') >= 0:
            type = 'ocean'
            if str.find('pacific') >= 0:
                name = 'pacific'
            elif str.find('georgia') >= 0:
                name = 'georgia'
            else:
                print 'Unknown name for open boundary %d:' % (ob+1)
                print str
                print '\nPlase check hgrid.gr3 and ensure all open boundaries '
                print 'have a name on the same line specifying the number of'
                print 'nodes for that boundary.'
                print 'e.g. "5 ! Nodes along boundary 4 ocean georgia'
                print 'Accepted ocean boundary names: georgia, pacific'
                name = '?'
        elif str.find('river') >= 0:
            type = 'river'
            if str.find('beaver') >= 0:
                name = 'beaver'
            elif str.find('bonneville') >= 0:
                name = 'bonneville'
            elif str.find('willamette') >= 0:
                name = 'willamette'
            elif str.find('morrison') >= 0:
                name = 'morrison'
            elif str.find('b_power1') >= 0:
                name = 'b_power1'
            elif str.find('b_power2') >= 0:
                name = 'b_power2'
            elif str.find('b_spillway') >= 0:
                name = 'b_spillway'
            elif str.find('fraser') >= 0:
                name = 'fraser'
            else:
                print 'Unknown name for open boundary %d:' % (ob+1)
                print str
                print '\nPlase check hgrid.gr3 and ensure all open boundaries '
                print 'have a name on the same line specifying the number of'
                print 'nodes for that boundary.'
                print 'e.g. "5 ! Nodes along boundary 4 river fraser'
                print 'Accepted river names: fraser, bonneville, beaver, morrison,'
                print '                      b_power1, b_power2, b_spillway'
                name = '?'
        elif str.find('land') >= 0:
            type = 'land'
            name = 'land'
        elif str.find('island') >= 0:
            type = 'island'
            name = 'island'
        else:
            print 'Type of open boundary for boundary %d is not specified:' % (ob+1)
            print str
            print '\nPlease check hgrid.gr3 and ensure all boundaries are'
            print 'specified as ocean, river, land, or island on the line '       
            print 'specifying the number of nodes along the boundary.'
            print 'e.g. "23 ! Nodes along boundary 3 river beaver"'
            type = '?'
            name = '?'
        _tmp = np.zeros((nNodesThisOB))
        # Loop over all nodes for this open boundary
        for obNode in range(nNodesThisOB):
            _tmp[obNode] = int(f.readline().rstrip("\n").split()[0])
        setattr(obNodes, "ob_%d_nodes" % (ob+1), _tmp)
        setattr(obNodes, "ob_%d_type" % (ob+1), type)
        setattr(obNodes, "ob_%d_name" % (ob+1), name)
    return obNodes

def readHGrid(path, gridonly=False):
    """ 
        Reads in an hgrid.gr3 file and places info into an object 
        
        grid - Object with grid information
    """
    try:
        f = open(path,"r")

        # Read the header information   
        header = f.readline()       
        [nElems, nNodes] = f.readline().rstrip("\n").split()
        nElems = int(nElems); nNodes = int(nNodes)

        # Get all the node information
        nodes = np.zeros((nNodes,4))
        for node in range(nNodes):
            [n, x, y, z] = f.readline().rstrip("\n").split()
            nodes[node,0] = int(n)+1
            nodes[node,1] = float(x)
            nodes[node,2] = float(y)
            nodes[node,3] = float(z)

        # Get the connectivity table
        elems = np.zeros((nElems,4))
        for elem in range(nElems):
            [n, x, a, b, c] = f.readline().rstrip("\n").split()
            elems[elem,0] = int(n)+1
            elems[elem,1] = int(a)
            elems[elem,2] = int(b)
            elems[elem,3] = int(c)

        # Get boundary node information
        if gridonly == False:
            obNodes = _readBndryNodes(f) # Open boundaries
            lbNodes = _readBndryNodes(f) # Land boundaries
            ibNodes = _readBndryNodes(f) # Island boundaries

        # Stuff everything into a generic object 
        grid = Object()
        grid.nodes = nodes
        grid.elems = elems
        if gridonly == False:
            grid.obNodes = obNodes
            grid.lbNodes = lbNodes
            grid.ibNodes = ibNodes
        return grid

    except IOError:
        raise
