#ifndef SHKZ_GRIDVISUALIZER3_CLOUD_INTERFACE_H
#define SHKZ_GRIDVISUALIZER3_CLOUD_INTERFACE_H

#include <shiokaze/graphics/graphics_engine.h>
#include <shiokaze/core/recursive_configurable_module.h>
#include <shiokaze/array/array3.h>
#include "../../src/octreeliquid/macoctreegrid3.h" //added

SHKZ_BEGIN_NAMESPACE

namespace macotreeliquid3_namespace { struct grid3; }
//using namespace macotreeliquid3_namespace;
class gridvisualizer3_cloud_interface : public recursive_configurable_module {
public:
    DEFINE_MODULE(gridvisualizer3_cloud_interface,"Grid Visualizer 3D cloud","GridVisualizercloud","cloud grid visualizer module")
    /**
     \~english @brief Draw density for cloud grid.
     @param[in] g Graphics engine.
     @param[in] grid cloud grid.
     \~japanese @brief オクツリーグリッドの密度を描画する。
     @param[in] g グラフィックスエンジン。
     @param[in] grid オクツリーグリッド。
     */
    virtual void draw_qc(graphics_engine &g, macotreeliquid3_namespace::grid3 &grid) const = 0;
    virtual void draw_qv(graphics_engine &g, macotreeliquid3_namespace::grid3 &grid) const = 0;
    virtual void draw_qr(graphics_engine &g, macotreeliquid3_namespace::grid3 &grid) const = 0;
private:
	virtual void initialize( const shape3 &shape, double dx ) = 0;
	virtual void initialize( const configurable::environment_map &environment ) override {
		//
		assert(check_set(environment,{"shape","dx"}));
		initialize(
			get_env<shape3>(environment,"shape"),
			get_env<double>(environment,"dx")
		);
	}
};

using gridvisualizer3_cloud_ptr = std::unique_ptr<gridvisualizer3_cloud_interface>;
using gridvisualizer3_cloud_driver = recursive_configurable_driver<gridvisualizer3_cloud_interface>;
//
SHKZ_END_NAMESPACE

#endif