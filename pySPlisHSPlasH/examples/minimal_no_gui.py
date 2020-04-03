import pysplishsplash as sph


def main():
    base = sph.Exec.SimulatorBase()
    base.init(useGui=False)
    base.setValueFloat(base.STOP_AT, 10.0) # Important to have the dot to denote a float
    base.run()


if __name__ == "__main__":
    main()
