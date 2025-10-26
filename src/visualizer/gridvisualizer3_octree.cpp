#ifndef SHKZ_GRIDVISUALIZER3_OCTREE_H
#define SHKZ_GRIDVISUALIZER3_OCTREE_H

#include <shiokaze/visualizer/gridvisualizer3_octree_interface.h>
#include "../octreeliquid/macoctreegrid3.h"

SHKZ_BEGIN_NAMESPACE

class gridvisualizer3_octree : public gridvisualizer3_octree_interface {
public:
    virtual void draw_density_oc(graphics_engine &g, macotreeliquid3_namespace::grid3 &grid) const {
    g.point_size(3.0);
    g.begin(graphics_engine::MODE::POINTS);
    
    // オクツリーグリッドのセルを巡回して密度を描画
	        grid.iterate_active_cells([&]( const cell_id3 &cell_id, int tid ) {
	            double density_value = grid.density[cell_id.index];
	            vec3d pos = grid.get_cell_position(cell_id);
	            g.color4(1.0,1.0,1.0,density_value);
	            g.vertex3v(pos.v);
	        });
    
    g.end();
    g.point_size(1.0);
}
};

SHKZ_END_NAMESPACE

#endif