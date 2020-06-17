# Restrictions

- When modifying simulation parameters this is the recommended structure, as modification will only work after `base.initSimulation()` has been called.
```python
base.initSimulation()
sim = sph.Simulation.getCurrent()
sim.setValue...()
base.runSimulation()
base.cleanup()
```
- `setValue...()` and `getValue...()` functions cannot accept vectors as arguments yet

