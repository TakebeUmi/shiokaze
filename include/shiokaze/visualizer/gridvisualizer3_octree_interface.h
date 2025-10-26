#ifndef SHKZ_GRIDVISUALIZER3_OCTREE_INTERFACE_H
#define SHKZ_GRIDVISUALIZER3_OCTREE_INTERFACE_H

#include "gridvisualizer3_interface.h"
#include "../../src/octreeliquid/macoctreegrid3.h" //added

SHKZ_BEGIN_NAMESPACE

namespace macotreeliquid3_namespace { struct grid3; }

class gridvisualizer3_octree_interface : public gridvisualizer3_interface {
public:
    DEFINE_MODULE(gridvisualizer3_octree_interface,"Grid Visualizer 3D Octree","GridVisualizerOctree","Octree grid visualizer module")

    /**
     \~english @brief Draw density for octree grid.
     @param[in] g Graphics engine.
     @param[in] grid Octree grid.
     \~japanese @brief オクツリーグリッドの密度を描画する。
     @param[in] g グラフィックスエンジン。
     @param[in] grid オクツリーグリッド。
     */
    virtual void draw_density_oc(graphics_engine &g, macotreeliquid3_namespace::grid3 &grid) const = 0;
};

using gridvisualizer3_octree_ptr = std::unique_ptr<gridvisualizer3_octree_interface>;

SHKZ_END_NAMESPACE

#endif