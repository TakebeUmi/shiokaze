import mitsuba
from mitsuba.core import *
from mitsuba.render import SceneHandler
from mitsuba.render import RenderQueue, RenderJob
import multiprocessing
from struct import *
import sys
import os
import time
import math
import ConfigParser, json
#
# Load variables
config = ConfigParser.ConfigParser()
config.read('../mesh/common.ini')
#
# Start up the scheduling system with one worker per local core
scheduler = Scheduler.getInstance()
for i in range(0, multiprocessing.cpu_count()):
	scheduler.registerWorker(LocalWorker(i, 'wrk%i' % i))
scheduler.start()

end = int(sys.argv[1])
xml_name = sys.argv[2]
#
SampleCount = config.getint('Common','SampleCount')
TargetPos = json.loads(config.get('Common','TargetPos'))
OriginPos = json.loads(config.get('Common','OriginPos'))
#
img_path = '../'+xml_name+'_img'
if not os.path.exists(img_path):
	os.system('mkdir '+img_path)

for frame in range(0,end+1):
	#
	png_path = img_path+'/'+str(frame)+'_'+xml_name+'.png'
	while( not os.path.exists(png_path)):
		#
		# Get a reference to the thread's file resolver
		fileResolver = Thread.getThread().getFileResolver()
		#
		# Interval for compling movie
		interval = 10
		#
		# Path to mesh files
		mesh_file = '../mesh/'+str(frame)+'_mesh.serialized'
		if xml_name == 'transparent':
			mesh_file = '../mesh/'+str(frame)+'_mesh_enclosed.serialized'
		#
		while( not os.path.exists(mesh_file)):
			time.sleep(1)
		#
		print 'Starting frame', frame
		#
		# Scene parameters
		paramMap = StringMap()
		paramMap['target_x'] = str(TargetPos[0])
		paramMap['target_y'] = str(TargetPos[1])
		paramMap['target_z'] = str(TargetPos[2])
		paramMap['origin_x'] = str(OriginPos[0])
		paramMap['origin_y'] = str(OriginPos[1])
		paramMap['origin_z'] = str(OriginPos[2])
		paramMap['up'] = '0, 1, 0'
		paramMap['sample_count'] = str(SampleCount)
		#
		paramMap['mesh_filename'] = mesh_file
		paramMap['solid_filename'] = '../mesh/static_solids/levelset_solid.serialized'
		scene = SceneHandler.loadScene(fileResolver.resolve(xml_name+'.xml'),paramMap)
		#
		# Append moving solids if exists
		transforms_file = '../mesh/'+str(frame)+'_transforms.dat'
		if os.path.exists(transforms_file):
			print 'Loading moving solids...'
			pmgr = PluginManager.getInstance()
			f = open(transforms_file,'rb')
			(number,) = unpack('I',f.read(4))
			print 'Adding '+str(number)+' moving solids...'
			for i in range(0,number):
				(x,y,z) = unpack('fff',f.read(12))
				(r,rx,ry,rz) = unpack('ffff',f.read(16))
				scene.addChild(pmgr.create(
					{	'type' : 'serialized',
						'maxSmoothAngle' : 30.0,
						'filename' : '../mesh/moving_solids/polygon_'+str(i)+'.serialized',
						'toWorld' : Transform.translate(Vector(x,y,z)) * Transform.rotate(Vector(rx,ry,rz),r),
						'bsdf' : {
							'type' : 'diffuse',
							'reflectance' : Spectrum(0.6)
						}
					}))
		#
		# Create a queue for tracking render jobs
		queue = RenderQueue()
		exr_path = img_path+'/'+str(frame)+'_'+xml_name+'.exr'
		scene.setDestinationFile(exr_path)
		#
		# Create a render job and insert it into the queue
		job = RenderJob('myRenderJob', scene, queue)
		job.start()
		#
		# Wait for all jobs to finish and release resources
		queue.waitLeft(0)
		queue.join()
		#
		# Print some statistics about the rendering process
		print(Statistics.getInstance().getStats())
		#
		# Tonemapping
		os.system('mtsutil tonemap -g 3 '+exr_path)
		#
		# Sometimes compile movie
		if frame > 0 and frame % interval == 0:
			video_path = img_path+'/'+xml_name+'.mp4'
			os.system('rm -rf '+video_path)
			os.system('avconv -r 60 -i '+img_path+'/%d_'+xml_name+'.png -b:v 120000k -c:v libx264 -pix_fmt yuv420p -loglevel panic -filter:v lutyuv="y=gammaval(1.3)" '+video_path)

