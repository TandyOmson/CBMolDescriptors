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
from prog_ops import get_flags_as_dict, xtb_opt, xtb_extract

from collections.abc import Iterable

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

                    attr = Descriptor(attr,i)
                    df.at[mol_name,i] = attr

        return df
    
    def unpack_df_descriptors(self,df):
        """ Unpacks all multidimensional descriptors in a dataframe
        """
        for i in df.columns:
            df[i] = df[i].apply(lambda x: x.unpack_descriptor() if isinstance(x, Descriptor) else x)

        return df
    
class Descriptor:
    """ Class for storing a descriptor and all its elements, before entry into dataframe
    """
    def __init__(self, descriptor, desc_name):
        self.contents = descriptor
        self.desc_name = desc_name
        self.desc_shape = np.shape(descriptor)
        if self.desc_shape != ():
            self.contents = np.array(self.contents)

    def unpack_descriptor(self):
        """ Unpack descriptor into a dictionary containing all elements of the descriptor and labels
        """ 
        desc_dict = {}

        if isinstance(self.contents, Iterable):
            ndim = len(self.desc_shape)
            if ndim == 1:
                for count,i in enumerate(self.contents):
                    desc_dict[f"{self.desc_name}_{count}"] = i
            elif ndim == 2:
                for i in range(self.desc_shape[0]):
                    for j in range(self.desc_shape[1]):
                        desc_dict[f"{self.desc_name}_{i}{j}"] = self.contents[i][j]
            elif ndim == 3:
                for i in range(self.desc_shape[0]):
                    for j in range(self.desc_shape[1]):
                        for k in range(self.desc_shape[2]):
                            desc_dict[f"{self.desc_name}_{i}{j}{k}"] = self.contents[i][j][k]
            elif ndim == 4:
                for i in range(self.desc_shape[0]):
                    for j in range(self.desc_shape[1]):
                        for k in range(self.desc_shape[2]):
                            for l in range(self.desc_shape[3]):
                                desc_dict[f"{self.desc_name}_{i}{j}{k}{l}"] = self.contents[i][j][k][l]

        else:
            desc_dict[self.desc_name] = self.contents

        return desc_dict
            

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

    df_unpacked = calc.unpack_df_descriptors(df)
    df_unpacked.to_pickle("./../results/descriptors.pkl")
    print(df_unpacked)

if __name__ == "__main__":
    main()
