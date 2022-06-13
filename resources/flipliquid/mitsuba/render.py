#!/usr/bin/env python
#
#	resources/flip/wscript
#
#	This is part of Shiokaze, a research-oriented fluid solver for computer graphics.
#	Created by Ryoichi Ando <rand@nii.ac.jp> on Sep 13, 2017.
#
#	Permission is hereby granted, free of charge, to any person obtaining a copy of
#	this software and associated documentation files (the "Software"), to deal in
#	the Software without restriction, including without limitation the rights to use,
#	copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the
#	Software, and to permit persons to whom the Software is furnished to do so,
#	subject to the following conditions:
#
#	The above copyright notice and this permission notice shall be included in all copies
#	or substantial portions of the Software.
#
#	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
#	INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
#	PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
#	HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
#	CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
#	OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
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

# Start up the scheduling system with one worker per local core
scheduler = Scheduler.getInstance()
for i in range(0, multiprocessing.cpu_count()):
	scheduler.registerWorker(LocalWorker(i, 'wrk%i' % i))
scheduler.start()

end = int(sys.argv[1])
name = sys.argv[2]
#
# Load variables
config = ConfigParser.ConfigParser()
config.read('../mesh/common.ini')
#
SampleCount = config.getint('Common','SampleCount')
TargetPos = json.loads(config.get('Common','TargetPos'))
OriginPos = json.loads(config.get('Common','OriginPos'))
LiquidColor = json.loads(config.get('Common','LiquidColor'))
#
img_path = '../'+name+'_img'
if not os.path.exists(img_path):
	os.system('mkdir '+img_path)
#
for frame in range(0,end+1):
	#
	png_path = img_path+'/'+str(frame)+'_'+name+'.png'
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
		if name == 'transparent':
			mesh_file = '../mesh/'+str(frame)+'_mesh_enclosed.serialized'

		particle_file = '../mesh/'+str(frame)+'_particles.dat'
		#
		while( not os.path.exists(mesh_file)):
			time.sleep(1)
		#
		while( not os.path.exists(particle_file)):
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
		paramMap['liquid_color'] = str(LiquidColor[0])+', '+str(LiquidColor[1])+', '+str(LiquidColor[2])
		#
		paramMap['mesh_filename'] = mesh_file
		paramMap['solid_filename'] = '../mesh/static_solids/levelset_solid.serialized'
		scene = SceneHandler.loadScene(fileResolver.resolve(name+'.xml'),paramMap)
		#
		# Append ballistic particles if exists
		if os.path.exists(particle_file):
			print 'Loading particles...'
			pmgr = PluginManager.getInstance()
			f = open(particle_file,'rb')
			(number,) = unpack('I',f.read(4))
			print 'Adding '+str(number)+' particles...'
			for i in range(0,number):
				(x,y,z,r,) = unpack('ffff',f.read(16))
				ParticleColor = Spectrum()
				ParticleColor.fromSRGB(LiquidColor[0],LiquidColor[1],LiquidColor[2])
				if name == 'mesh':
					scene.addChild(pmgr.create(
									   {'type' : 'sphere',
									   'center' : Point(x,y,z),
									   'radius' : r,
									   'ref' : {
									   'type' : 'plastic',
									   'diffuseReflectance' : ParticleColor,
									   'intIOR' : 1.5
									   }
									   }))
				elif name == 'transparent':
					scene.addChild(pmgr.create(
									   {'type' : 'sphere',
									   'center' : Point(x,y,z),
									   'radius' : r,
									   'ref' : {
									   'type' : 'dielectric',
									   'intIOR' : 1.5
									   }
									   }))
		#
		# Create a queue for tracking render jobs
		queue = RenderQueue()
		scene.setDestinationFile(img_path+'/'+str(frame)+'_'+name+'.exr')
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
					{	'type' : 'ply',#'serialized',
						#'maxSmoothAngle' : 30.0,
						'filename' : '../mesh/moving_solids/polygon_'+str(i)+'.ply',#serialized',
						'toWorld' : Transform.translate(Vector(x,y,z)) * Transform.rotate(Vector(rx,ry,rz),r),
						#'bsdf' : {
						#	'type' : 'diffuse',
						#	'reflectance' : Spectrum(0.6)
						#}
					}))
		#
		# Create a render job and insert it into the queue
		job = RenderJob('myRenderJob',scene,queue)
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
		exr_path = img_path+'/'+str(frame)+'_'+name+'.exr'
		os.system('mtsutil tonemap -g 3 '+exr_path)
		#
		# Sometimes compile movie
		if frame > 0 and frame % interval == 0:
			video_path = img_path+'/'+name+'.mp4'
			os.system('rm -rf '+video_path)
			os.system('avconv -r 60 -i '+img_path+'/%d_'+name+'.png -b:v 120000k -c:v libx264 -pix_fmt yuv420p -loglevel panic -filter:v lutyuv="y=gammaval(1.3)" '+video_path)