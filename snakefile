


OUTMAIN = config["outmain"]
PRENAME=config["prename"]

SAMPLES=list(config["samples"].keys())
BAMS=[config["samples"][x] for x in SAMPLES]

SATC = config["satc"]
PARFINDER=config["parfinder"]

include: "rules/input_prepping.smk"
include: "rules/satc.smk"

