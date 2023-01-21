-- Lua (Keep this comment, this is an indication for editor's 'run' command)

S = scene_graph.create_object({classname='OGF::MeshGrob',name='S'})
S.I.Shapes.create_cube({x1=0,y1=0,z1=0,x2=1,y2=1,z2=1})
S.I.Shapes.create_sphere({precision=4,radius=0.5,center='0 0 0'})
S.I.Shapes.create_sphere({precision=4,radius=0.5,center='0 0 1'})
S.I.Shapes.create_sphere({precision=4,radius=0.5,center='0 1 0'})
S.I.Shapes.create_sphere({precision=4,radius=0.5,center='0 1 1'})
S.I.Shapes.create_sphere({precision=4,radius=0.5,center='1 0 0'})
S.I.Shapes.create_sphere({precision=4,radius=0.5,center='1 0 1'})
S.I.Shapes.create_sphere({precision=4,radius=0.5,center='1 1 0'})
S.I.Shapes.create_sphere({precision=4,radius=0.5,center='1 1 1'})

S.I.Surface.triangulate()
S.I.Experiment.intersect_surface()

