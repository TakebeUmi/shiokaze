
#include <shiokaze/visualizer/gridvisualizer3_cloud_interface.h>
#include "../octreeliquid/macoctreegrid3.h"
#include <shiokaze/core/console.h>
#include <vector>
SHKZ_USING_NAMESPACE
using namespace macoctreeliquid3_namespace;
class gridvisualizer3_cloud : public gridvisualizer3_cloud_interface {
protected:
    virtual void draw_qc(graphics_engine &g, macotreeliquid3_namespace::grid3 &grid) const {
    g.point_size(3.0);
    g.begin(graphics_engine::MODE::POINTS);
	auto func = [&]( const macotreeliquid3_namespace::cell_id3 &cell_id ) {
		double qc_value = grid.qc[cell_id.index];
		vec3d pos = grid.get_cell_position(cell_id);
		g.color4(1.0,0.5,0.5,qc_value*100.0);
		g.vertex3v(pos.v); 
		//if (qc_value > 0.000001) console::dump("%u:%f at (%f,%f,%f)\n", (unsigned)cell_id.index, qc_value,pos[0],pos[1],pos[2]);
	};
	
    //std::vector<vec3d> points;
    // オクツリーグリッドのセルを巡回して密度を描画
				// grid.iterate_active_cells([&]( const macotreeliquid3_namespace::cell_id3 &cell_id, int tid ) {
				// 	double qc_value = grid.qc[cell_id.index];
				// 	vec3d pos = grid.get_cell_position(cell_id);
				// 	//points.push_back(pos);
				// 	g.color4(1.0,1.0,1.0,qc_value);
				// 	vec3d pos_ = {0.1,0.1,0.1};
				// 	g.vertex3v(pos_.v); 
				// 	console::dump("%u:%f at (%f,%f,%f)\n", (unsigned)cell_id.index, qc_value,pos[0],pos[1],pos[2]);
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

	virtual void draw_qv(graphics_engine &g, macotreeliquid3_namespace::grid3 &grid) const {
    g.point_size(3.0);
    g.begin(graphics_engine::MODE::POINTS);
	auto func = [&]( const macotreeliquid3_namespace::cell_id3 &cell_id ) {
		double qv_value = grid.qv[cell_id.index];
		vec3d pos = grid.get_cell_position(cell_id);
		g.color4(0.5,1.0,0.5,qv_value*10);
		g.vertex3v(pos.v); 
		//if (qc_value > 0.000001) console::dump("%u:%f at (%f,%f,%f)\n", (unsigned)cell_id.index, qc_value,pos[0],pos[1],pos[2]);
	};
	
    //std::vector<vec3d> points;
    // オクツリーグリッドのセルを巡回して密度を描画
				// grid.iterate_active_cells([&]( const macotreeliquid3_namespace::cell_id3 &cell_id, int tid ) {
				// 	double qc_value = grid.qc[cell_id.index];
				// 	vec3d pos = grid.get_cell_position(cell_id);
				// 	//points.push_back(pos);
				// 	g.color4(1.0,1.0,1.0,qc_value);
				// 	vec3d pos_ = {0.1,0.1,0.1};
				// 	g.vertex3v(pos_.v); 
				// 	console::dump("%u:%f at (%f,%f,%f)\n", (unsigned)cell_id.index, qc_value,pos[0],pos[1],pos[2]);
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

    virtual void draw_qr(graphics_engine &g, macotreeliquid3_namespace::grid3 &grid) const {
    g.point_size(3.0);
    g.begin(graphics_engine::MODE::POINTS);
	auto func = [&]( const macotreeliquid3_namespace::cell_id3 &cell_id ) {
		double qr_value = grid.qr[cell_id.index];
		vec3d pos = grid.get_cell_position(cell_id);
		g.color4(0.5,0.5,1.0,qr_value);
		g.vertex3v(pos.v); 
		//if (qc_value > 0.000001) console::dump("%u:%f at (%f,%f,%f)\n", (unsigned)cell_id.index, qc_value,pos[0],pos[1],pos[2]);
	};
	
    //std::vector<vec3d> points;
    // オクツリーグリッドのセルを巡回して密度を描画
				// grid.iterate_active_cells([&]( const macotreeliquid3_namespace::cell_id3 &cell_id, int tid ) {
				// 	double qc_value = grid.qc[cell_id.index];
				// 	vec3d pos = grid.get_cell_position(cell_id);
				// 	//points.push_back(pos);
				// 	g.color4(1.0,1.0,1.0,qc_value);
				// 	vec3d pos_ = {0.1,0.1,0.1};
				// 	g.vertex3v(pos_.v); 
				// 	console::dump("%u:%f at (%f,%f,%f)\n", (unsigned)cell_id.index, qc_value,pos[0],pos[1],pos[2]);
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
		config.get_bool("DrawqcOC",m_param.draw_qc,"Should draw qc OC");
		config.get_bool("DrawqvOC",m_param.draw_qv,"Should draw qv OC");
		config.get_bool("DrawqrOC",m_param.draw_qr,"Should draw qr OC");
	}

	shape3 m_shape;
	double m_dx;
	//
	struct Parameters {
		bool draw_qc {true};	
		bool draw_qv {true};
		bool draw_qr {true};
	};
	//
	Parameters m_param;
};

extern "C" module * create_instance() {
	return new gridvisualizer3_cloud();
}
//
extern "C" const char *license() {
	return "MIT";
}