#ifndef H_OGF_MESH_TOOLS_MESH_GROB_NEIGH_TOOLS_H
#define H_OGF_MESH_TOOLS_MESH_GROB_NEIGH_TOOLS_H

#include <OGF/Experiment/common/common.h>
#include <OGF/mesh_gfx/tools/mesh_grob_tool.h>

namespace OGF {

    gom_attribute(help, "select facet neighborhoods")
    gom_attribute(message, "btn1: select; btn3: unselect")
    gom_class Experiment_API MeshGrobNeighborTool : public MeshGrobTool {
    public:
        MeshGrobNeighborTool(ToolsManager* parent);
        void grab(const RayPick& p_ndc) override;
        void reset() override;
    };

}

#endif
