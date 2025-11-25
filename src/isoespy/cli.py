# handler function
# src/isoespy/cli.py
import sys

from .isoespy_de import main as de_main
from .isoespy_diu import main as diu_main
from .isoespy_ff import main as ff_main
from .isoespy_edger import main as edger_main
from .run_DRIMSeq_stageR import main as drimseq_main
from .isoespy_makefa import main as makefa_main
from .isoespy_makegtf import main as makegtf_main
from .isoespy_orfpred import main as orfpred_main


COMMANDS = {
    "de": de_main,
    "diu": diu_main,
    "ff": ff_main,
    "edger": edger_main,
    "drim": drimseq_main,
    "fas": makefa_main,
    "gtf": makegtf_main,
    "orf": orfpred_main,
}


def _print_help() -> None:
    print("Usage: isoespy <command> [options]\n")
    print("Available commands:")
    print("$ isoespy  de [options]       Differential transcript expression visualization")
    print("isoespy de visualizes the exon structures of each isoform of the target gene along with the expression distributions of the two groups and the results of the differential expression (DE) analysis.")
    print("")

    print("$ isoespy  diu [options]      Differential isoform usage visualization")
    print("isoespy diu visualizes the exon structures of each isoform of the target gene, the isoform usage proportions within each group, and the changes in usage between the two groups.")
    print("")

    print("$ isoespy  ff [options]       Functional feature visualization")
    print("isoespy ff visualizes the exon structures and functional features of each isoform of the target gene.")
    print("")

    print("$ isoespy  edger [options]    edgeR-based transcript-level DE analysis")
    print("isoespy edger is a wrapper for edgeR. This command performs transcript-level differential expression analysis using the edgeR R package.")
    print("")

    print("$ isoespy  drim [options]     DRIMSeq and stageR-based DIU analysis")
    print("isoespy drim is a wrapper for DRIMSeq and stageR. It performs gene-level differential isoform usage (DIU) analysis using these R packages.")
    print("")

    print("$ isoespy  fas [options]      Handle FASTA file")
    print("isoespy fas extracts FASTA sequences from a GTF File. To extract nucleotide sequences for transcript models from an input GTF file, specify feature exon and type nucleotide. If you wish to obtain amino acid sequences instead, specify feature CDS and type amino_acid. Because sequence retrieval is performed using the pysam library, the reference genome index file (.fai) must be located in the same directory as the reference FASTA file.")
    print("")

    print("$ isoespy  gtf [options]      Handle GTF file")
    print("isoespy gtf integrates GTFprep annotations into a GTF file. By providing a GTFprep file that contains functional annotation information, isoespy gtf generates a new GTF file in which these annotations are appended to the corresponding transcripts.")
    print("")

    print("$ isoespy  orf [options]      ORF prediction")
    print("isoespy orf predicts ORFs for transcripts that lack CDS annotations in GTF file. For any transcript in the input GTF that does not contain a CDS feature, isoespy uses CPAT and TransDecoder to predict coding sequences and appends the corresponding CDS annotations to the GTF. The hexamer and model options in the following command correspond to human. When analyzing other organisms, please provide the appropriate pretrained datasets for the target species.")
    print("")

def main() -> int:
    """Entry point for the `isoespy` command with subcommands."""

    # サブコマンドが指定されていない / help のとき
    if len(sys.argv) < 2 or sys.argv[1] in ("-h", "--help"):
        _print_help()
        return 0

    cmd = sys.argv[1]

    # サブコマンド名を argv から取り除き、
    # 残りをそのまま各 main() に渡せるようにする
    sys.argv = [sys.argv[0]] + sys.argv[2:]

    if cmd in COMMANDS:
        # 各 isoespy_*.py の main() は従来通り sys.argv から引数を読むはずなので、
        # ここでは単純に呼び出すだけでよい
        return COMMANDS[cmd]()
    else:
        print(f"Unknown command: {cmd}\n")
        _print_help()
        return 1


if __name__ == "__main__":
    raise SystemExit(main())


