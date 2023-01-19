-- Lua (Keep this comment, this is an indication for editor's 'run' command)

S = scene_graph.create_object({classname='OGF::MeshGrob',name='S'})
S.I.Shapes.create_cube({x1=0,y1=0,z1=0,x2=1,y2=1,z2=1})
S.I.Shapes.create_cube({x1=1,y1=0.3,z1=0.3,x2=2,y2=0.7,z2=0.7})
S.I.Surface.triangulate()
S.I.Experiment.intersect_surface()

