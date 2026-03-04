

rule find_par:
    input:
        tsv = expand(os.path.join(OUTMAIN, "satc_wins", "output", PRENAME + "_win{winsize}", PRENAME + "_winClassification.tsv"),
                    winsize=config["winsizes"])


rule do_satc:
    input:
        scaffList = os.path.join(OUTMAIN, "satc", PRENAME+"_scaffSex.tsv"),
        xzlist = os.path.join(OUTMAIN, "satc", PRENAME+"_XZ_scaff.list")


rule run_satc:
    input:
        idx_list = os.path.join(OUTMAIN, "idxstats", PRENAME+"_all_idxstats.list")
    output:
        scaffList = os.path.join(OUTMAIN, "satc", PRENAME+"_scaffSex.tsv"),
        xzlist = os.path.join(OUTMAIN, "satc", PRENAME+"_XZ_scaff.list"),
        samplesex =  os.path.join(OUTMAIN, "satc", PRENAME+"_sampleSex.tsv")
    params:
        outprefix = os.path.join(OUTMAIN, "satc", PRENAME),
        minlen = config["minlen"]
    conda:
        "../envs/satc.yml"
    resources:
        mem_mb = 10000,
        runtime = 300
    shell: """
        "{conda_env}/bin/Rscript" {SATC} -i {input.idx_list} -o {params.outprefix} --minLength {params.minlen}
    """



rule run_satc_findPar:
    input:
        depmat = os.path.join(OUTMAIN, "satc_wins", "inputs", PRENAME + "_goodauto_XZ_scaff_win{winsize}.depths.txt"),
        samplesex =  os.path.join(OUTMAIN, "satc", PRENAME + "_sampleSex.tsv")
    output:
        pngfolder = directory(os.path.join(OUTMAIN, "satc_wins", "output", PRENAME + "_win{winsize}", PRENAME + "_plots")),
        tsv = os.path.join(OUTMAIN, "satc_wins", "output", PRENAME + "_win{winsize}", PRENAME + "_winClassification.tsv")
    params:
        outprefix = os.path.join(OUTMAIN, "satc_wins", "output", PRENAME + "_win{winsize}", PRENAME)
    conda:
        "../envs/satc.yml"
    resources:
        mem_mb = 10000,
        runtime = 300
    shell: """
        mkdir -p {output.pngfolder}
        "{conda_env}/bin/Rscript" {PARFINDER} {input.depmat} {input.samplesex} {params.outprefix}
    """
