""" Main for molecular descriptor calculation program 

    Flags are set to calculate descriptors in config.ini
    Sections in the config file must correspond to the name of the module, format <sectionname>_desc.py

    Usage:
        python main.py <inputfile>
"""

import configparser
import importlib

import pandas as pd
import numpy as np

import sys
import os
import pprint

from inout import get_csv_smiles, get_sdfs
from prog_ops import get_flags_as_dict, xtb_calc, xtb_opt, xtb_extract, df_item

class configManager:
    """ Class for managing config.ini file
    
        Attributes:
            config (configparser): config.ini file
            flags (dict): flags for calculating descriptors
    """

    def __init__(self, config_file):
        self.config = configparser.ConfigParser(allow_no_value=True)
        self.config.read(config_file)
        self.section_names = self.config.sections()

    # FLAGS
    def get_section_flags(self, section):
        """ Get flags for a section in config.ini """
        runflag = self.config.getboolean(section, "run")
        if runflag == True:
            flags = get_flags_as_dict(self.config,section)
        else:
            flags = {}

        return runflag, flags

    # I/O
    def get_input(self):
        """ Get input file path """
        pass

    def get_output(self):
        """ Get output file path """
        pass

    def get_print_options(self):
        """ Get print options """
        pass

class DescriptorCalculator:
    """ Class for calculating molecular descriptors
        Imports descriptor calculation modules based on flags
        Each descriptor function should take an rdkit mol object, and return a descriptor value
        The calc_desc class methods then imports the relevant methods form the modules and calls the descriptor function on the dataframe
    """

    def __init__(self, configfile):
        self.config = configManager(configfile)
        # Get section names
        sections = self.config.section_names
        self.runflags = {}
        for i in sections:
            runflag, flagdict = self.config.get_section_flags(i)
            setattr(self, f"run_{i}", runflag)
            setattr(self, f"flags_{i}", flagdict)
            self.runflags[i] = runflag

    # Calc section batchwise (add piecewise method?)
    def calc_section(self, df, section):
        module = importlib.import_module(f"{section}_desc")
        section_flags = getattr(self, f"flags_{section}")
        
        for mol, mol_name in zip(df["guestmol"],df.index):
            for i,j in section_flags.items():
                if j and i != "run":
                    try:
                        attr = getattr(module, i)(mol)
                    except:
                        attr = np.nan
                    if not isinstance(attr, list):
                        df.at[mol_name,i] = attr
                    elif isinstance(attr, list) or isinstance(attr, np.array):
                        # Add to dataframe as class object for ease of reading
                        attr_df = df_item(attr)
                        df.at[mol_name,i] = attr_df

        return df

def main():
    """ Main function for calculating molecular descriptors """
    inputpath = sys.argv[1]

    # Get input file type
    if inputpath.endswith('.smi'):
        df = get_csv_smiles(inputpath)
    else:
        df = get_sdfs(inputpath)

    print(df)

    pp = pprint.PrettyPrinter(indent=4)

    # Initialize DescriptorCalculator
    calc = DescriptorCalculator("config.ini")
    print("\nSections to be run:\n==============")
    pp.pprint(calc.runflags)
    # Run all sections whose section flag is true
    for section,runflag in calc.runflags.items():
        if runflag == True:
            if section == "xtb":
                # Run an xtb calculation, add the dictionary of results as a property of the guestmol object
                if os.path.isdir("./xtb_temp"):
                    df = xtb_extract(df)
                else:
                    df = xtb_opt(df)
                print(df)
            df = calc.calc_section(df,section)

    print(df)
    df.to_pickle("/home/spine/DProjects/DDescriptCalc/results/descriptors.pkl")

if __name__ == "__main__":
    main()