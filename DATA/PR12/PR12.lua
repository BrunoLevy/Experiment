-- Lua  

scene_graph.clear()
S = scene_graph.create_object('OGF::MeshGrob','S')
E = S.I.Editor

function R(a,x,y,z)
    x = x - 0.5
    y = y - 0.5
  s = math.sin(a)
  c = math.cos(a)
  return s*x+c*y, -c*x+s*y, z
end

function cube(a)
 
  o = E.nb_vertices
  
  E.create_vertex(R(a, 0.0, 0.0, 0.0))
  E.create_vertex(R(a, 0.0, 0.0, 1.0))
  E.create_vertex(R(a, 0.0, 1.0, 0.0))
  E.create_vertex(R(a, 0.0, 1.0, 1.0))
  E.create_vertex(R(a, 1.0, 0.0, 0.0))
  E.create_vertex(R(a, 1.0, 0.0, 1.0))
  E.create_vertex(R(a, 1.0, 1.0, 0.0))
  E.create_vertex(R(a, 1.0, 1.0, 1.0))

  E.create_quad(o+2,o+3,o+1,o+0)
  E.create_quad(o+4,o+5,o+7,o+6)
  E.create_quad(o+0,o+1,o+5,o+4)
  E.create_quad(o+1,o+3,o+7,o+5)
  E.create_quad(o+3,o+2,o+6,o+7)
  E.create_quad(o+2,o+0,o+4,o+6)
end 

for i = 0,2 do
   cube(i*0.1)
end

--S.I.Experiment.intersect_surface({check_constraints=true})
--S.save('result.geogram')


