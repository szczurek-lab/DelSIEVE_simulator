# DelSIEVE_simulator
The simulator for DelSIEVE is modified based on [SIEVE_simulator](https://github.com/szczurek-lab/SIEVE_simulator). It is used in the DelSIEVE paper to generate simulated datasets.

## Compilation

DelSIEVE_simulator uses cmake to build the executable.

```bash
$ mkdir build
$ cd build
$ cmake ..
$ cmake --build .
```

The built executable is located: `build/bin/DelSIEVE_simulator`.

## Run

To run DelSIEVE_simulator, use the example configuration file `parameters`:

```bash
$ ./build/bin/SIEVE_simulator -Fparameters
```

## Configuration files

DelSIEVE_simulator uses a similar input as CellCoal does. Please refer to the manuals of CellCoal. Those configuration files used in the DelSIEVE paper can be found in the [DelSIEVE_benchmark_pipeline](https://github.com/szczurek-lab/DelSIEVE_benchmark_pipeline) repository, under directory `simulation_configs`.

