#include <OGF/Experiment/tools/mesh_grob_neighbor_tools.h>

namespace OGF {
    
    MeshGrobNeighborTool::MeshGrobNeighborTool(
        ToolsManager* parent
    ) : MeshGrobTool(parent) {
    }
    
    void MeshGrobNeighborTool::reset() {
        Object_var I = mesh_grob()->query_interface("Selections");
        if(!I.is_null()) {
            I->invoke_method("show_facets_selection");
        }
    }
    
    void MeshGrobNeighborTool::grab(const RayPick& p_ndc) {
        Attribute<bool> selection(
            mesh_grob()->facets.attributes(), "selection"
        );
        index_t f = pick_facet(p_ndc);
        if(f != index_t(-1)) {
            bool on = (p_ndc.button == 1);
            selection[f] = on;
            Logger::out("Neigh") << "Picked " << f << std::endl;
            for(
                index_t le = 0;
                le < mesh_grob()->facets.nb_vertices(f); ++le
            ) {
                index_t f_adj = mesh_grob()->facets.adjacent(f,le);
                if(f_adj != index_t(-1)) {
                    selection[f_adj] = on;
                    Logger::out("Neigh") << " neighbor " << f_adj << std::endl;
                }
            }
        }
        mesh_grob()->update();
    }
}

