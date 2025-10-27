
#include <shiokaze/visualizer/gridvisualizer3_octree_interface.h>
#include "../octreeliquid/macoctreegrid3.h"
#include <shiokaze/core/console.h>
#include <vector>
SHKZ_USING_NAMESPACE
using namespace macoctreeliquid3_namespace;
class gridvisualizer3_octree : public gridvisualizer3_octree_interface {
protected:
    virtual void draw_density_oc(graphics_engine &g, macotreeliquid3_namespace::grid3 &grid) const {
    g.point_size(3.0);
    g.begin(graphics_engine::MODE::POINTS);
	auto func = [&]( const macotreeliquid3_namespace::cell_id3 &cell_id ) {
		double density_value = grid.density[cell_id.index];
		vec3d pos = grid.get_cell_position(cell_id);
		g.color4(1.0,1.0,1.0,density_value);
		g.vertex3v(pos.v); 
		console::dump("%u:%f at (%f,%f,%f)\n", (unsigned)cell_id.index, density_value,pos[0],pos[1],pos[2]);
	};
    //std::vector<vec3d> points;
    // オクツリーグリッドのセルを巡回して密度を描画
				// grid.iterate_active_cells([&]( const macotreeliquid3_namespace::cell_id3 &cell_id, int tid ) {
				// 	double density_value = grid.density[cell_id.index];
				// 	vec3d pos = grid.get_cell_position(cell_id);
				// 	//points.push_back(pos);
				// 	g.color4(1.0,1.0,1.0,density_value);
				// 	vec3d pos_ = {0.1,0.1,0.1};
				// 	g.vertex3v(pos_.v); 
				// 	console::dump("%u:%f at (%f,%f,%f)\n", (unsigned)cell_id.index, density_value,pos[0],pos[1],pos[2]);
				// });
		
	for( char depth=0; depth<grid.layers.size(); ++depth ) {
		const auto &layer = *grid.layers[depth];
		layer.active_cells.const_serial_actives([&]( int i, int j, int k, const auto &it ) {
			func({depth,vec3i(i,j,k),it()});
		});
	}
	// vec3d pos_ = {0.1,0.1,0.1};
	// g.color4(1.0,1.0,1.0,1.0);
	// g.vertex3v(pos_.v);
	// g.color4(1.0,0.0,0.0,1.0);
	// pos_ = {0.2,0.2,0.2};
	// g.vertex3v(pos_.v);
    g.end();
    g.point_size(1.0);
	}

	virtual void initialize( const shape3 &shape, double dx ) override {
		m_shape = shape;
		m_dx = dx;
	}
	virtual void initialize( const filestream &file ) override {
		file.r(m_shape);
		file.r(m_dx);
	}
	virtual void serialize( const filestream &file ) const override {
		file.w(m_shape);
		file.w(m_dx);
	}
	virtual void configure( configuration &config ) override {
		config.get_bool("DrawDensityOC",m_param.draw_density_oc,"Should draw density OC");
	}

	shape3 m_shape;
	double m_dx;
	//
	struct Parameters {
		bool draw_density_oc {true};		
	};
	//
	Parameters m_param;
};

extern "C" module * create_instance() {
	return new gridvisualizer3_octree();
}
//
extern "C" const char *license() {
	return "MIT";
}