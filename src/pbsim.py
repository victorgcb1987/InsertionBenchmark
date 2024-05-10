from subprocess import run


def run_pbsim(strategy="wgs", depth=100, min_length=1000, max_length=30000,
              method="errhmm", method_model=None, reference=None, sequencing_preset="ccs"):
    cmd = ["pbsim"]
    if strategy in ["wgs"]:
        cmd.append("--strategy {}".format(strategy))
    else:
        raise ValueError("{} is not a valid strategy")
    cmd.append("--depth {}".format(str(depth)))
    cmd.append("--length-min {}".format(str(min_length)))
    cmd.append("--length-max {}".format(str(max_length)))
    if method in ["errhmm", "qshmm"]:
        if method_model.exists():
            cmd.append("--method {} --{} {}".format(method, method, str(method_model)))
        else:
            raise RuntimeError("Model file not found:".format(str(method_model)))
    else:
        raise ValueError("{} is not a valid method".format(method))
    if reference.exists():
        cmd.append("--genome {}".format(str(reference)))
    else:
        raise RuntimeError(("Sequence file not found: {}".format(str(reference))))
    if sequencing_preset != "ccs":
        prefix = "sd"
        run(" ".join(cmd), shell=True)
    else:
        prefix = "ccs"
        cmd.append("--prefix ccs --pass-num 10  --seed 100")
        run(" ".join(cmd), shell=True)
        run("samtools view -bS -@ 40 ccs_0001.sam > ccs_0001.bam", shell=True)
        run("ccs all -j 40 ccs_0001.bam ccs_0001.fastq", shell=True)
    run("gzip -c {}_0001.fastq > {}_0001.fastq.gz".format(prefix, prefix), shell=True)
    run("rm {}_0001.fastq".format(prefix), shell=True)
    run("gzip -c {}_0001.maf > {}_0001.maf.gz".format(prefix), shell=True)
    run("rm {}_0001.maf".format(prefix), shell=True)


    #~/soft/pbsim3/src/pbsim --strategy wgs --depth 100  --length-min 10000  --length-max 30000 --method qshmm --qshmm ~/soft/pbsim3/data/QSHMM-RSII.model --depth 100 --prefix ccs --pass-num 10  --seed 100  --genome  Nuclear_with_insertions.fasta