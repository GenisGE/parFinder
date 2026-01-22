
import sys
import yaml
import glob
import os

species = sys.argv[1]
bamfolder=sys.argv[2]
# species = "Sturnus_vulgaris"
prename = species.split("_")[0].lower()

outfile = "config_" + prename + ".yaml"

if os.path.isfile(outfile):
    print("config file already exists, will exit")
    raise SystemExit


#bamfolder="/path/to/folder/with/bam/files/"
##if organized by sps cna use this
#mainfolder = "/path/to/folder/with/species/bam/folders/"
#bamfolder = mainfolder + species + "/BAM/" 



bams = glob.glob(bamfolder+"*.qc.bam")
sample_names = [os.path.basename(x).split(".")[0] for x in bams]

if len(bams) == 0:
    print("no bams in folder for "+ species + ", will exit without config")
    raise SystemExit


mydict = {
    "satc": "SATC/satc.R",
    "parfinder": "scripts/parFinder.R",
    "outmain": "results/" + species,
    "prename": prename,
    "samples": dict(zip(sample_names, bams)),
    "winsizes": [50000, 100000],
    "minlen": 100000
} 


fh = open(outfile, "w+")
yaml.dump(mydict, fh)

print("created configfile " + outfile)