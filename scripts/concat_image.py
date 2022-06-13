import os, sys

if len(sys.argv) > 2:
		index = 0
		concat_label = ""
		for name in sys.argv[1:]:
			concat_label += name
		while True:
				command = "convert"
				exit_loop = False
			for name in sys.argv[1:]:
					filename = name + "_" + str(index) + ".png"
					if not os.path.exists(filename):
						exit_loop = True
						break
				command += " " + filename
				if exit_loop:
					break
				command += " -append " + concat_label + "_" + str(index) + ".png"
				print command
				os.system(command)
				index += 1

