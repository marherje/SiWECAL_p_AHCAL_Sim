# Masking

The `masker` utility to apply the masking from data to the simulation. The file containing the masked channels has the usual format (layer chip chn_list). Examples can be found in the eos at `/eos/project-s/siw-ecal/TB2022-06/beamData/semi_online_mon_data/raw_siwecal_[RUN_NUMBER]/masked_channels.txt`.


## Example

```bash
root -l
root [0] .L masker.C
root [1] masker("raw_siwecal_90474_masked_channels.txt", "/data_ilc/flc/jimenez/simulations/TB2022-06/CONF6/build/ECAL_QGSP_BERT_conf6_e-_8GeV_5kevt_build.root")
```