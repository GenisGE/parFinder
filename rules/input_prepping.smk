

rule do_depmat:
    input:
        outmat = expand(os.path.join(OUTMAIN, "satc_wins", "inputs", PRENAME+ "_goodauto_XZ_scaff_win{winsize}.depths.txt"),
                        winsize = config["winsizes"])



rule do_idxstats:
    input:
        bam = lambda wildcards: config["samples"][wildcards.s]
    output:
        idxstats = os.path.join(OUTMAIN, "idxstats", "{s}.idxstats")
    conda:
        "../envs/samtools.yml"
    resources:
        mem_mb = 50000,
        runtime = 300
    threads: 3
    shell:
        "samtools idxstats {input.bam} > {output.idxstats}"


rule do_idxstats_list:
    input:
        expand(os.path.join(OUTMAIN, "idxstats", "{s}.idxstats"),
                s = config["samples"].keys())
    output:
        idxlist = os.path.join(OUTMAIN, "idxstats", PRENAME+"_all_idxstats.list")
    resources:
        mem_mb = 1000,
        runtime = 100
    run:
        with open(output.idxlist, "w+") as fh:
            for idxstat in input:
                fh.write(idxstat + "\n")


rule make_windows_scafflist:
    # select the largest well-behaved autosomal (more strict, requrie p-value of difference between male and female > 0.1) for normalization of windows,
    # plus all "XZ" scaffolds to find the par region in
    input:
        scaff_sex = os.path.join(OUTMAIN, "satc", PRENAME+"_scaffSex.tsv"),
        xz_list = os.path.join(OUTMAIN, "satc", PRENAME+"_XZ_scaff.list")
    output:
        scaff_chromsize = os.path.join(OUTMAIN, "satc_wins", "inputs", PRENAME+"_goodauto_XZ_scaff.chromsize.txt")
    resources:
        mem_mb = 2000,
        runtime = 100
    shell: """
    echo "will start running"
    set +o pipefail
    awk '$3 == "FALSE" && $4 == "FALSE" && $5 > 0.1' {input.scaff_sex} | sort -hrk2 | head -1 | cut -f1,2 > {output.scaff_chromsize}
    awk '$3 == "TRUE" || $4 == "TRUE"' {input.scaff_sex} | sort -hrk2  | cut -f1,2 >> {output.scaff_chromsize}
    """


rule make_windows_bed:
    input:
        scaff_chromsize = os.path.join(OUTMAIN, "satc_wins", "inputs", PRENAME+ "_goodauto_XZ_scaff.chromsize.txt")
    output:
        bedwins = os.path.join(OUTMAIN, "satc_wins", "inputs", PRENAME+ "_goodauto_XZ_scaff_win{winsize}.bed")
    resources:
        mem_mb = 5000,
        runtime = 1000
    conda:
        "../envs/bedtools.yml"
    shell: 
        """
        bedtools makewindows -g {input.scaff_chromsize} -w {wildcards.winsize} > {output.bedwins}
        """


rule make_window_depmat:
    input:
        bams = BAMS,
        bedwins = os.path.join(OUTMAIN, "satc_wins", "inputs", PRENAME+ "_goodauto_XZ_scaff_win{winsize}.bed")
    output:
        outmat = os.path.join(OUTMAIN, "satc_wins", "inputs", PRENAME+ "_goodauto_XZ_scaff_win{winsize}.depths.txt")
    params:
        header = "scaffold start end len " + " ".join(SAMPLES)
    resources:
        mem_mb = 50000,
        runtime = 2000
    conda:
        "../envs/samtools.yml"
    shell: """
        echo {params.header} > {output.outmat}
        samtools bedcov {input.bedwins} {input.bams} | awk '{{$3=($3" "($3-$2)); print $0}}'|  tr -s " " >> {output.outmat}
    """
