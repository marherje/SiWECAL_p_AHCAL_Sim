#!/bin/env python3

## Example usage:
## ./prepare_Shaping.py --template ../../steering/templates/Shaping_TB2022-03.xml
##      --output_filename ../../steering/Shaping/TB2022-03/ECAL_QGSP_BERT_conf1_e-_3GeV_30.xml
##      --LCIOInputFiles /data_ilc/flc/jimenez/simulations/TB2022-03/CONF1/lcio/ECAL_QGSP_BERT_conf1_e-_3GeV_30.slcio
##      --LCIOOutputFile /data_ilc/flc/jimenez/simulations/TB2022-03/CONF1/lcio/ECAL_QGSP_BERT_conf1_e-_3GeV_30_shaped.slcio
##      --MaxRecordNumber 1000

import xmltodict
import sys, os, argparse, pathlib
from pprint import pprint as pp 


MAX_REC = 1000

parser = argparse.ArgumentParser(description="Write Shaping processor steering files from template")
parser.add_argument('--template', help="Template to start from")
parser.add_argument('--output_filename', help="Output xml filename")
parser.add_argument('--LCIOInputFiles', nargs='+', help="Input LCIO file(s)")
# parser.add_argument('--LCIOOutputFile', nargs='+', help="Output LCIO file")
parser.add_argument('--LCIOOutputFile', help="Output LCIO file")
parser.add_argument('--MaxRecordNumber', type=int, default=MAX_REC, help="N events (default 1k)")
#TODO should be a list
parser.add_argument('--Input_Collections', default="SiEcalCollection", help="Input collection to use (default: SiEcalCollection)")
parser.add_argument('--Output_Collections', default="ShapedSiEcalCollection", help="Output collection to use (default: ShapedSiEcalCollection)")
parser.add_argument('--tuples', help="String of tuples, in the format \"((key_path_1, val1), (key_path_2, val2), ...)\" to overwrite keys in the steering (not implemented yet)")

def transformer(args, f):
    d = xmltodict.parse(f.read())
    # pp(d)
    # Should be done in a better way, but it was annoying me
    for ipar, par in enumerate(d["marlin"]["global"]["parameter"]):
        if par["@name"] == "LCIOInputFiles":
            d["marlin"]["global"]["parameter"][ipar]["#text"] = f'\n\t\t'+ f'\n\t\t'.join(args.LCIOInputFiles) + f'\n\t\t'
        elif par["@name"] == "MaxRecordNumber":
            d["marlin"]["global"]["parameter"][ipar]["@value"] = args.MaxRecordNumber
    
    for iproc, proc in enumerate(d["marlin"]["processor"]):
        if proc["@name"] == "MyLCIOOutputProcessor":
            for ipar, par in enumerate(proc["parameter"]):
                if par["@name"] == "LCIOOutputFile":
                    d["marlin"]["processor"][iproc]["parameter"][ipar]["#text"] = args.LCIOOutputFile
        if proc["@name"] == "SiWECALShaping":
            for ipar, par in enumerate(proc["parameter"]):
                if par["@name"] == "Input_Collections":
                    d["marlin"]["processor"][iproc]["parameter"][ipar]["#text"] = args.Input_Collections
                elif par["@name"] == "Output_Collections":
                    d["marlin"]["processor"][iproc]["parameter"][ipar]["#text"] = args.Output_Collections

    
    
    # #for key, val in (("LCIOInputFiles", args.LCIOInputFiles),
    # #                 ("MaxRecordNumber", str(args.MaxRecordNumber))):
    # for ipar, param_key in enumerate(f_dict["marlin"]["global"]["parameter"]):
    #     if param_key["@name"] == "LCIOInputFiles":
    #         f_dict["marlin"]["global"]["parameter"][ipar]["#text"] = f'\n\t\t'+ f'\n\t\t'.join(args.LCIOInputFiles) + f'\n\t\t'
    #     elif param_key["@name"] == "LCIOOutputFile":
    #         #f_dict["marlin"]["global"]["parameter"][ipar]["#text"] = f'\n\t\t'+ f'\n\t\t'.join(args.LCIOOutputFile) + f'\n\t\t'
    #         f_dict["marlin"]["global"]["parameter"][ipar]["#text"] = args.LCIOOutputFile
    #     elif param_key["@name"] == "MaxRecordNumber":
    #         f_dict["marlin"]["global"]["parameter"][ipar]["@value"] = str(args.MaxRecordNumber)
    
    # for iproc, proc in enumerate(f_dict["marlin"]["processor"]):
    #     print("This proc", proc, '\n\n')
    #     if proc["@name"] == "MyLCIOOutputProcessor":
    #         for ipar, param_key in enumerate(proc["parameter"]):
    #             if param_key["@name"] == "LCIOOutputFile":
    #                 f_dict["marlin"]["processor"]["parameter"][ipar]["#text"] = args.Output_Collections
    # #TODO Implement tuples here
    # #print(f_dict)


    return(xmltodict.unparse(d, pretty=True))

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

