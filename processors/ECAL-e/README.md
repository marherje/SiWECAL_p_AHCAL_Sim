# Processors for digitization

## Notes on the processors

The main processors of interest for simulation and digitization are:

- `ConversionProcessor`: To make Raw energy -> MIP conversion.
- `ShapingProcessor`: Simulates the pulse shaping in the ASICs (in progress of adaptation for TB21/22).
- `ShapedConversion`: A new conversion to MIP after the effect of the shaping (in progress of adaptation for TB21/22).

There are other processors for misc tasks such as `(Digi)LCIO2BuildProcessor` for transforming LCIO (Sim)CalorimeterHits to the "build" format used in data.

## Setup and running

```
source /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/init_ilcsoft.sh # Or your favorite ilcsoft version
./build.sh
export MARLIN_DLL=$PWD/lib/libDigitization.so:$MARLIN_DLL
Marlin <path_to_steering_file.xml>
```

Under `run_scripts` there are per-processor `prepare_[processor].py` (python3) scripts to help with the steering file writing, to ease multiple file processing. Input templates can be found under `steering/templates/`. Caveat: these scripts need to be run in an environment with python3 (the source of ilcsoft above is not py3, so, you should use a different one).

### Examples

Suppose we have lcio files for simulated mu- (MIP, for conversion) and e- for TB2022-06. Then, prepare the conversion xml file (py3 env):

```bash
./run_scripts/prepare_Conversion.py --template steering/templates/Conversion.xml --output_filename steering/Conversion/TB2022-06_CONF6_mu-_conv_test.xml --LCIOInputFiles [mu-_LCIO_FILES] --ConvAuxFile [OUTPUT_CONV_ROOT_FILENAME] --MaxRecordNumber 5000 --tb_conf TB2022-06_CONF6 --MIPFitMode 1
```

Run the processor with (ilcsoft env):

```bash
Marlin steering/Conversion/TB2022-06_CONF6_mu-_conv_test.xml
```

The `[OUTPUT_CONV_ROOT_FILENAME]` is a file that contains the resulting histograms and fit for the conversion per layer. The processor itself reports the conversion factors, but you can use the following script to get the conversion factors that should be added to processors steering file templates downstream (py3):

```bash
./run_scripts/conversion2steering.py [OUTPUT_CONV_ROOT_FILENAME]
```

For deriving a ROOT file out of the SLCIO (with the conversion), prepare the steering file with (py3):

```bash
 ./run_scripts/prepare_LCIO2Build.py --template steering/templates/LCIO2Build.xml
      --output_filename steering/LCIO2Build/TB2022-06_CONF6_e-_10GeV_test.xml
      --LCIOInputFiles [e-_LCIO_FILES]
      --OutputBuildFile [OUTPUT_DERIVED_ROOT_FILENAME]
      --MaxRecordNumber 5000
```

Finally, run the processor as usual (ilcsoft env):

```bash
Marlin steering/Conversion/TB2022-06_CONF6_e-_10GeV_test.xml
```

REMEMBER of CONF6:

"TB2022-06_CONF6" : [ [4.2, True, 650],
                          [4.2, True, 650],
                          [4.2, True, 650],
                          [4.2, True, 650],
                          [4.2, True, 500],
                          [4.2, True, 500],
                          [4.2, True, 500],
                          [4.2, True, 500],
                          [5.6, True, 500],
                          [5.6, True, 500],
                          [5.6, True, 320],
                          [5.6, True, 320],
                          [5.6, True, 320],
                          [5.6, True, 320],
                          [5.6, True, 320], # Careful, wf1 0.5
                        ],