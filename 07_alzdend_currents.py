"""
Modified protocol for a Purkinje cell subject to atrophied dendrites
in Alzheimer's. This protocol can be called with command line arguments
as follows:
python 07_alzdend_currents.py cutoff protocol
# Examples:
python 07_alzdend_currents.py 250	# call the function to use a cell with max dendritic length of 250 micron

python 07_alzdend_currents.py 250 1 # call the function and specify to connect several current clamps to dendrites
									# instead of one to soma. It will maintain a consistent density of current clamps
									# so that full size cells have almost 500 and smaller cells have proportionally fewer

# There is the option to run this interactively, with the -i flag:
python -i 07_alzdend_currents.py

# And then, if you pause and restart the simulation in the runcontrol panel, after you are finished with it you can
# call (again) the function that writes results. It will overwrite the previous results with the longer duration data:

>>> saveresults(cell,cutlength)
"""

from smallPurkinje import Purkinje, smPurk
from neuron import h
import multiprocessing
import numpy as np
import sys, math

#fixed time step only
Fixed_step = h.CVode()
Fixed_step.active(0) #the model does not work with the variable time step!

#Instantiation of the cell template
compcell = Purkinje()

stimdata = dict()
stimdata['stim0del'] = 100
stimdata['stim0dur'] = 4000
stimdata['stim0amp'] = 2 #0.01 #2

stimdata['timeglobal'] = 4000

# EXTRA STUFF
TM = 0
maxDlength = compcell.getmaxlength() # may be able to get this programmatically in NEURON
#cutlength=maxDlength - (1-TM)*(130 - 85)*maxDlength/130
cutlength = 250
if len(sys.argv)>1:
	cutlength   = float(sys.argv[1]) # can pass in cutlength as an argument at the command line, ex:
									 # python -i 07_alzdend_currents.py 250

print("cutlength = ", cutlength)
compcell.prune(cutoff=cutlength)

cell = smPurk()
compcell = None
h.PlotShape().printfile("PictureOfCell_"+str(cutlength)+".eps")
h.PlotShape().show(0)
h.PlotShape().printfile("PictureOfCell_"+str(cutlength)+"_diam.eps")


#this code discover the number of cores available in a CPU and activate the multisplit to use them all.
cores = multiprocessing.cpu_count()
h.load_file("parcom.hoc")
p = h.ParallelComputeTool()
p.change_nthread(cores,1)
p.multisplit(1)
print('cores', cores)

#Neuron control menu
h.nrncontrolmenu()

#Voltage graph
h('load_file("vm.ses")')

stimtype=0
if len(sys.argv)>2:
	stimtype=1

if stimtype==0:
	#ORIGINAL CODE
	stim = [h.IClamp(0.5,sec=cell.soma)]
else:
	#NEW STUFF
	import random

	random.seed(1)

	l = []
	stim = []
	test = 0

	L=cell.alldendlength()
	print('combined length of all dendrites remaining = ', round(L))
	
	# num syns to add (constant density along dendrites)
	numSyn=math.ceil(L*0.04)
	
	for x in range(0, numSyn): # place that number of synapses on dendrites
		test = random.randint(0, len(cell.dend)) # randomly pick which dendrites to place on
		stim.append(h.IClamp(0.5,sec=cell.dend[test])) # place in middle of dendrite
		l.append(test) # keep track of where placed
		
		# apply same current injection at each place:
		stim[x].delay = stimdata['stim0del']
		stim[x].dur = stimdata['stim0dur']
		stim[x].amp = stimdata['stim0amp']*0.0001
		
	#for i in l:
	#	print(i)
	print('place ', numSyn, ' inputs')
	##END NEW STUFF

#Basic properties of the simulation. dt, temperature, sim duration and initial voltage
h.dt = 0.025
h.celsius = 37
h.tstop = stimdata['timeglobal']
h.v_init = -65

#initialization and run.    
def initialize():
	h.finitialize()
	h.run()
	h.PlotShape().show(0)

def saveresults(cell,cutlength): # you can call this again manually if you run more time in the runcontrol
	#save vectors for time and voltage
	time = np.array(cell.rec_t)
	if len(time)==0:
		time=np.arange(0,h.t+h.dt,h.dt)
	vm_soma = np.array(cell.vm_soma)
	vm_NOR3 = np.array(cell.vm_NOR3)
	if len(time)>len(vm_soma):
		time=np.arange(0,h.t,h.dt)

	timevolt_soma = np.column_stack((time,vm_soma))
	timevolt_NOR3 = np.column_stack((time,vm_NOR3))

	#save files
	np.savetxt('07_vm_soma_' + str(cutlength) + '_stim' + str(stimtype) + '.txt', timevolt_soma, delimiter = ' ')
	np.savetxt('07_vm_NOR3_' + str(cutlength) + '_stim' + str(stimtype) + '.txt', timevolt_NOR3, delimiter = ' ')
	print('saved results in: ' + '07_vm_soma_' + str(cutlength) + '_stim' + str(stimtype) + '.txt')

	
initialize()
saveresults(cell,cutlength)