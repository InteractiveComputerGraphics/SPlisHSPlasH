# Embedded Python

## Build with embedded Python support

To enable the embedded Python support just activate the CMake option USE_EMBEDDED_PYTHON which is by default turned off. Please ensure that CMake finds the Python interpreter. This can be achieved by setting the PYTHON_EXECUTABLE to the file path of the python interpreter.

## Run simulator with embedded Python support

Make sure that the environment variables PYTHONHOME and PYTHONPATH are set to the directory of your Python installation. Also make sure that the pythonXX.dll (where XX defines the version) is in your path.

When running the simulator with embedded Python support, the new tab 'Scripts' should appear in the GUI. Here you can load a script file or reload it. 

Alternatively, you can load the Python script directly in a scene. Just use the json key ```scriptFile``` in the scene file to define the location of the script. If a relative path is used, the simulator assumes that it is relative to the scene file. 

## Writing a script

First, you have to import the module ```splishsplash```, e.g. by: ```import splishsplash as sph ```

When the script is loaded, the simulator will call the function ```init()``` automatically. If you defined a function ```step()```, it will be called in each simulation step. If you defined a function ```reset()```, it will be called when the simulation is reset. Moreover, you can define additional functions that can be called using the GUI.  

### init(base)

If this function is defined, it is called automatically when the script is loaded or reloaded. The parameter base contains the current SimulationBase object.

### step()

If this function is defined, it is called automatically in each simulation step. 

### reset()

If this function is defined, it is called automatically in each simulation reset. 

### Additional commands

Additional commands that can be executed via the GUI (buttons will be added) can be defined by a list of strings. The list must be called ```function_list``` and must contain the names of functions that should be called, e.g.
```python
function_list = ['command', 'command2']
```

## Example 

This is a simple example script which prints the current simulation time, the position of particle 0 and a counter value in each step:

```python
import splishsplash as sph
import numpy as np

counter = 0
function_list = ['command', 'command2']

def init(base):
    global counter
    print("init test")
    counter = 1
    
def step():
    global counter
    sim = sph.Simulation.getCurrent()
    fluid = sim.getFluidModel(0)
    tm = sph.TimeManager.getCurrent()
    print(fluid.getPosition(0))
    print(tm.getTime())
    print(counter)
    counter += 1
    print("---")
    
def reset():
    print("reset test")
    
def command():
    print("tst cmd")
        
def command2():
    print("tst cmd2")
```