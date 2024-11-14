# bed-pileup.py v0.1

# Owen Chapman
# 03/09/2022
# Command line usage:
# python3 bed-pileup.py --dir [directory-containing-bed-files] -o [output-bed-file]

# WARNING: Do not combine bed files originating from different genome assemblies!
# the results will be meaningless.

import pandas as pd
import pyranges as pr
import os
import argparse
import sys

class BedDict(dict):
    path=None
    def __init__(self,path):
        self.path=path
        for entry in  os.listdir(self.path):
            f = os.path.join(self.path,entry)
            if os.path.isfile(f) and f.endswith('.bed'):
                try:
                    self[entry]=pr.read_bed(f)
                except:
                    print(f"could not read file: {str(f)}")

    def dict_pileup(self):
        beds = pr.count_overlaps(self)
        beds = beds.as_df()
        beds["Name"] = beds[beds.columns.to_list()[3:]].sum(axis=1)
        beds = beds[["Chromosome","Start","End","Name"]]
        beds = beds[beds.Name != 0] # why are there zero entries at all??
        return pr.PyRanges(beds)

if __name__ == "__main__":
    # Parse input params
    parser = argparse.ArgumentParser('Generate .bedgraph by stacking bed files.')
    parser.add_argument('-d','--dir',help='Directory containing .bed files')
    parser.add_argument('-o','--outfile', help='outfile (optional)',default=None)
    args = parser.parse_args()

    if args.outfile == None:
        #raise NotImplementedError()
        write=lambda x: sys.stdout.write(x)
        print('Writing to stdout')
    else:
        write=lambda x: x
        print('Writing to',args.outfile)

    # Run
    beds = BedDict(args.dir)
    bed = beds.dict_pileup()

    # Write
    bed = bed.as_df()
    write(bed.to_csv(args.outfile,sep='\t',header=False,index=False))



