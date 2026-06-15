import json
from importlib.metadata import version

from .runCommands import command


def getVersionsAutodatabase(args):
    """ get python package and software versions for afanc-autodatabase
    """

    version_dict = {
        "afanc" : version("afanc"),
        "Biopython" : version("Biopython"),
        "numpy" : version("numpy"),
        "pysam" : version("pysam"),
        "pandas" : version("pandas"),
        "scipy" : version("scipy"),
        "ncbi_taxonomy_date" : args.ncbi_date
    }

    version_dict["mash"] = _get_command_version(["mash", "--version"])

    version_dict["Kraken2"] = _parse_kraken_version(_get_command_version(["kraken2", "-v"]))

    versions_json = f"{args.autoDB_WDir}/versions.json"

    with open(versions_json, "w") as fout:
        json.dump({ "afanc-autodatabase_versions" : version_dict }, fout, indent = 4)


def getVersionsScreen():
    """ Gets python package and softeware versions for afanc-screen
    """

    version_dict = {
        "afanc" : version("afanc"),
        "Biopython" : version("Biopython"),
        "numpy" : version("numpy"),
        "pysam" : version("pysam"),
        "pandas" : version("pandas"),
        "scipy" : version("scipy")
    }

    version_dict["Kraken2"] = _parse_kraken_version(_get_command_version(["kraken2", "-v"]))

    version_dict["BWA"] = _get_command_version("bwa 2>&1 | head -n 1", first_nonempty_line=True)
    version_dict["samtools"] = _get_command_version(["samtools", "--version"], first_nonempty_line=True)
    version_dict["samclip"] = _get_command_version(["samclip", "--version"], first_nonempty_line=True)
    version_dict["freebayes"] = _get_command_version(["freebayes", "--version"], first_nonempty_line=True)
    version_dict["bcftools"] = _get_command_version(["bcftools", "--version"], first_nonempty_line=True)

    return version_dict


def _get_command_version(argv, first_nonempty_line=False):
    try:
        stdout, stderr = command(argv, "GET-VERSIONS").run_comm_quiet(1)
    except Exception as exc:
        return f"unavailable: {exc}"

    text = stdout.decode(errors="replace").strip()
    if not text:
        text = stderr.decode(errors="replace").strip()

    if first_nonempty_line:
        for line in text.splitlines():
            if line.strip():
                return line.strip()
        return ""

    return text


def _parse_kraken_version(version_text):
    if "version " in version_text:
        return version_text.split("\n")[0].split("version ", 1)[1]
    return version_text
