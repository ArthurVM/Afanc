from collections import defaultdict
from importlib.metadata import version

from Afanc.utilities.runCommands import command

def get_versions_screen():
    """ Gets python package versions for afanc-screen
    """

    version_dict = {
        "afanc" : version("afanc"),
        "Biopython" : version("Biopython"),
        "numpy" : version("numpy"),
        "pysam" : version("pysam"),
        "pandas" : version("pandas"),
        "scipy" : version("scipy")
    }

    mash_stdout, mash_stderr = command(f"mash --version", "GET-VERSIONS").run_comm_quiet(1)
    version_dict["mash"] = mash_stdout.decode().strip("\n")

    k2_stdout, k2_stderr = command(f"kraken2 -v", "GET-VERSIONS").run_comm_quiet(1)
    version_dict["Kraken2"] = k2_stdout.decode().split("\n")[0].split("version ")[1]

    bt2_stdout, bt2_stderr = command(f"bowtie2 --version", "GET-VERSIONS").run_comm_quiet(1)
    version_dict["Bowtie2"] = bt2_stdout.decode().split("\n")[0].split("version ")[1]

    return version_dict
