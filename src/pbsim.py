from subprocess import run


def run_pbsim(strategy="wgs", depth=100, min_length=1000, max_length=30000,
              method="errhmm", method_model=None, reference=None):
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
    run(" ".join(cmd), shell=True, capture_output=True)