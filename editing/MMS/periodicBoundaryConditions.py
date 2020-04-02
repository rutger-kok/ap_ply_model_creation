# Periodic Boundary Conditions

vps = session.viewports.values()[0]
modelName = vps.displayedObject.modelName

activeModel = mdb.models[modelName]
specimenPart = activeModel.parts['Specimen']
specimenEdges = specimenPart.edges

xmax = ymax = 37.5
xmin = ymin = -12.5
xmid = (xmax+xmin)/2.0 + xmin
ymid = (ymax+ymin)/2.0 + ymin
tol = 0.001
# Creating nodes set at the edges
leftEdge = specimenEdges.getByBoundingBox(xmin-tol, ymin-tol, -1.0, xmin+tol,
                                          ymax+tol, 1.0)
leftEdgeSet = specimenPart.Set(edges=leftEdge, name='Left Edge')

rightEdge = specimenEdges.getByBoundingBox(xmax-tol, ymin-tol, -1.0, xmax+tol,
                                           ymax+tol, 1.0)
rightEdgeSet = specimenPart.Set(edges=rightEdge, name='Right Edge')

# Creating individual nodes
# Nodes left edge
a = []
for i in leftEdgeSet.nodes:
    a = a+[(i.coordinates[1], i.label)]
a.sort()
for ind, j in enumerate(a):
    leftNodeSet = specimenPart.Set(name='LeftNode-'+str(ind),
                                   nodes=specimenPart.nodes[(j[1]-1):(j[1])])

# Nodes right edge
b = []
for k in rightEdgeSet.nodes:
    b = b+[(k.coordinates[1], k.label)]
b.sort()
for ind2, m in enumerate(b):
    rightNodeSet = specimenPart.Set(name='RightNode-'+str(ind2),
                                    nodes=specimenPart.nodes[(m[1]-1):(m[1])])


# Applying constraint equation to each node for PBC
for num in range(1, len(leftNodeSet)+1):
    leftnodename = 'Specimen Instance.LeftNode-{}'.format(num)
    rightnodename = 'Specimen Instance.RightNode-{}'.format(num)
    activeModel.Equation(name='Constraint-{}'.format(num),
                         terms=((1.0, leftnodename, 1),
                                (-1.0, rightnodename, 1)))
