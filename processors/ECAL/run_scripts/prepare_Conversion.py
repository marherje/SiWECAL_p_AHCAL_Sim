#!/bin/env python3

## Example usage:
## ./prepare_Conversion.py --template ../../steering/templates/Conversion.xml
##      --output_filename TB2022-06_CONF6_mu-_150GeV_conversion_test.xml
##      --LCIOInputFiles /data_ilc/flc/jimenez/simulations/TB2022-08/CONF6/lcio/ECAL_QGSP_BERT_conf6_mu-_150GeV_{10..19}.slcio
##      --ConvAuxFile /data_ilc/flc/jimenez/simulations/TB2022-06/CONF6/build/ECAL_QGSP_BERT_conf6_mu-_150GeV_5kevt_conversion.root
##      --tb_conf TB2022-06_CONF6
##      --MipFitMode 1
##      --MaxRecordNumber 5000

import xmltodict
import sys, os, argparse, pathlib
from confs import W_Confs


MAX_REC = 5000

parser = argparse.ArgumentParser(description="Write processor steering files from a template")
parser.add_argument('--template', help="Template to start from")
parser.add_argument('--output_filename', help="Output xml filename")
parser.add_argument('--LCIOInputFiles', nargs='+', help="Input LCIO file(s)")
parser.add_argument('--MaxRecordNumber', type=int, default=MAX_REC, help="N events (default 5k)")
#TODO should be a list
parser.add_argument('--Input_Collections', default="SiEcalCollection", help="Input collection to use (default: SiEcalCollection)")
parser.add_argument('--ConvAuxFile', help="Output auxiliary root file for conversion")
parser.add_argument('--tb_conf', help="TB configuration name in confs.py")
parser.add_argument('--MIPFitMode', help="Fit mode: 1 (Gaus), 2 (Landau), 3 (LanGaus)")
parser.add_argument('--tuples', help="String of tuples, in the format \"((key_path_1, val1), (key_path_2, val2), ...)\" to overwrite keys in the steering (not implemented yet)")

def transformer(args, f):
    f_dict = xmltodict.parse(f.read())
    #for key, val in (("LCIOInputFiles", args.LCIOInputFiles),
    #                 ("MaxRecordNumber", str(args.MaxRecordNumber))):
    for ipar, param_key in enumerate(f_dict["marlin"]["global"]["parameter"]):
        if param_key["@name"] == "LCIOInputFiles":
            #f_dict["marlin"]["global"]["parameter"][ipar]["#text"] = args.LCIOInputFiles 
            f_dict["marlin"]["global"]["parameter"][ipar]["#text"] = f'\n\t\t'+ f'\n\t\t'.join(args.LCIOInputFiles) + f'\n\t\t'
        elif param_key["@name"] == "MaxRecordNumber":
            f_dict["marlin"]["global"]["parameter"][ipar]["@value"] = str(args.MaxRecordNumber)

    # print(W_Confs[args.tb_conf])
    if args.tb_conf in W_Confs: thicknesses = " ".join([str(conf[2]) for conf in W_Confs[args.tb_conf]])
    for key, val in (("Input_Collections", args.Input_Collections),
                     ("ConvAuxFile", args.ConvAuxFile),
                     ("MIPFitMode", args.MIPFitMode),
                     ("SiThicknesses", thicknesses)): 
        for ipar, param_key in enumerate(f_dict["marlin"]["processor"]["parameter"]):
            if param_key["@name"] == key:
                f_dict["marlin"]["processor"]["parameter"][ipar]["#text"] = val 
    #TODO Implement tuples here
    #print(f_dict)
    return(xmltodict.unparse(f_dict, pretty=True))

# argparse

if __name__ == "__main__":
    args = parser.parse_args()
    filename = args.template
    
    # Roundtrip xml - json - xml
    with open(filename) as f1:
        new_xml = transformer(args, f1)
    new_filename = args.output_filename
    with open(new_filename, "w") as f2:
        f2.write(new_xml)
    print("New xml file written in ", new_filename)
    print("Run with:\nMarlin", new_filename)