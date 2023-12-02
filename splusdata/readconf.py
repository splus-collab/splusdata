import argparse
import splusdata
import pandas as pd
import os

from yaml import load, dump
try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper
    
operations = [
    'stamp',
    'lupton_rgb',
    'trilogy_image'
    #'field_frame',
    #'checkcoords',
    #'query'
]

def handle_operation_type(operation):
    if isinstance(operation["type"], str):
        ops = [operation["type"]]
    for op in ops:
        if op not in operations:
            raise ValueError("Operation type not supported: {}".format(op))
    if "file_path" in operation:
        f = pd.read_csv(operation["file_path"])

        if "ra" not in f.columns:
            raise ValueError("File {} does not have column 'ra'".format(operation["file_path"]))
        if "dec" not in f.columns:
            raise ValueError("File {} does not have column 'dec'".format(operation["file_path"]))
        
        
        
    return ops


def main():
    parser = argparse.ArgumentParser(description='splusdata - Download SPLUS catalogs, FITS and more')
    parser.add_argument('config_file', metavar='config_file', type=str, help='Configuration file')
    
    print(os.getcwd())
    
    args = parser.parse_args()

    configfile = os.path.join(os.getcwd(), args.config_file)

    stream = open(configfile, 'rb')
    
    data = load(stream, Loader=Loader)
    print(data)
    print(type(data))

    try:
        conn = splusdata.Core(data["user"], data["password"])
    except:
        conn = splusdata.Core()

        

if '__name__' == '__main__':
    main()  