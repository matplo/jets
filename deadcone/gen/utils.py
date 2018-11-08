import os

def grab_module_environment():
	rets = []
	module_home = os.environ.get('MODULESHOME')
	init_scipt = os.path.join(module_home, 'init/bash')
	if os.path.isfile(init_scipt):
		rets.append('. {}'.format(init_scipt))
	module_path = os.environ.get('MODULEPATH')
	splits = module_path.split(':')
	for sp in splits:
		rets.append('module use {}'.format(sp))
	module_loaded = os.environ.get('LOADEDMODULES')
	splits = module_loaded.split(':')
	for sp in splits:
		rets.append('module load {}'.format(sp))
	return '\n'.join(rets)
