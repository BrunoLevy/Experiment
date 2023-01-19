-- PR1 crashes
S=scene_graph.create_object({classname='OGF::Mesh',name='S'})
S.I.Shapes.create_icosahedron()
S.I.Shapes.create_sphere()
S.I.Experiment.intersect_surface()
