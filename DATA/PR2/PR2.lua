-- Lua PR2 
scene_graph.clear()
S=scene_graph.create_object({classname='OGF::Mesh',name='S'})
S.I.Shapes.create_icosahedron()
S.I.Shapes.create_sphere({precision=2})
-- S.I.Experiment.intersect_surface()

