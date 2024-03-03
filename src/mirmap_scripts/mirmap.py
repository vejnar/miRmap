#!/usr/bin/env python3

#
# Copyright © 2024 Charles E. Vejnar
#
# This is free software, licensed under the GNU General Public License v3.
# See /LICENSE for more information.
#

import argparse
import itertools
import os
import shutil
import sys

import mirmap.if_lib_spatt
import mirmap.if_lib_viennarna
import mirmap.scores
import mirmap.target


table_columns_annots = (
    "mirna_id",
    "transcript_stable_id",
    "target_id",
    "seed_length",
    "mirna_start",
    "mirna_end",
    "seed_start",
    "seed_end",
)
table_columns_scores = (
    "tgs_au",
    "tgs_position",
    "tgs_pairing3p",
    "tgs_score",
    "dg_duplex",
    "dg_binding",
    "dg_duplex_seed",
    "dg_binding_seed",
    "dg_open",
    "dg_total",
    "prob_exact",
    "prob_binomial",
    "cons_bls",
    "selec_phylop",
    "mirmap_score",
)
table_columns_annots_1to1 = (
    "mirna_id",
    "transcript_stable_id",
    "num_target",
    "num_seed6",
    "num_seed7",
)


def check_exe(names):
    """Check executable is found."""
    for name in names:
        if shutil.which(name) is None:
            raise FileNotFoundError(f"{name} missing")


def read_seq(seq=None, fname_fasta=None, fname_tab=None, id=None):
    """Read sequence(s) from FASTA or TSV file."""
    assert seq is not None or fname_fasta is not None or fname_tab is not None, "FASTA or TSV input required"
    if seq is not None:
        if id is not None:
            return {id: seq}
        else:
            return {"user": seq}
    else:
        if fname_fasta is not None:
            seqs = mirmap.utils.load_fasta(fname_fasta)
        elif fname_tab:
            with open(fname_tab, "rt") as f:
                seqs = dict(line.rstrip().split("\t") for line in f)
        if id is not None:
            return {k: v for k, v in seqs.items() if k == id}
        else:
            return seqs


def format_number(v):
    """Format number controlling the number of decimals reported."""
    if isinstance(v, int):
        return str(v)
    else:
        s1 = str(v)
        s2 = f"{v:.8}"
        if len(s1) < len(s2):
            return s1
        else:
            return s2


def main(argv=None):
    """Main."""
    # Parameters
    if argv is None:
        argv = sys.argv
    parser = argparse.ArgumentParser(description="Predict miRNA targets.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-m", "--mirna", dest="mirna_seq", action="store", help="miRNA sequence")
    parser.add_argument("-n", "--mirna-id", dest="mirna_id", action="store", help="miRNA IDs")
    group.add_argument("-a", "--mirna-fasta", dest="mirna_fasta", action="store", help="miRNA Fasta file")
    group.add_argument("-b", "--mirna-tab", dest="mirna_tab", action="store", help="miRNA tabulated file")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-t", "--transcript", dest="transcript_seq", action="store", help="Transcript sequence")
    parser.add_argument("-i", "--transcript-id", dest="transcript_id", action="store", help="Transcript IDs")
    group.add_argument(
        "-f", "--transcript-fasta", dest="transcript_fasta", action="store", help="Transcript Fasta file"
    )
    group.add_argument(
        "-u", "--transcript-tab", dest="transcript_tab", action="store", help="Transcript tabulated file"
    )
    parser.add_argument(
        "-g",
        "--aggregate",
        dest="aggregate",
        action="store_true",
        help="Aggregate multiple targets (miRNA-mRNA 1 to 1 relationships)",
    )
    parser.add_argument(
        "-s", "--aln", dest="path_aln", action="store", default=".", help="Path to multiple sequence alignment(s)"
    )
    parser.add_argument(
        "-d", "--mod", dest="path_mod", action="store", default=".", help="Path to evolutionary model(s)"
    )
    parser.add_argument("-e", "--tree", dest="path_tree", action="store", help="Path to Newick species tree")
    parser.add_argument("-o", "--output", dest="output", action="store", default="-", help="Output")
    parser.add_argument("-x", "--output-1to1", dest="output_1to1", action="store", default="-", help='Output "1 to 1"')
    parser.add_argument("-p", "--pretty-output", dest="pretty_output", action="store_true", help="Pretty output")
    parser.add_argument(
        "--path-libspatt2",
        dest="path_libspatt2",
        action="store",
        default="libspatt2.so",
        help="Path to libspatt2.so library",
    )
    parser.add_argument(
        "--path-phylofit",
        dest="path_phylofit",
        action="store",
        default="phyloFit",
        help="Path to the phyloFit executable",
    )
    parser.add_argument(
        "--path-phylop", dest="path_phylop", action="store", default="phyloP", help="Path to the phyloP executable"
    )
    parser.add_argument(
        "--temperature", dest="temperature", action="store", default=37.0, help="RNA folding temperature"
    )
    args = parser.parse_args(argv[1:])

    # Check
    exes = [args.path_libspatt2, args.path_phylop]
    if args.path_tree is not None:
        exes.append(args.path_phylofit)
    check_exe(exes)

    # Sequences input
    mirnas = read_seq(args.mirna_seq, args.mirna_fasta, args.mirna_tab, args.mirna_id)
    transcripts = read_seq(args.transcript_seq, args.transcript_fasta, args.transcript_tab, args.transcript_id)

    # Read tree
    if args.path_tree is not None:
        tree = open(args.path_tree, "rt").read()
    else:
        tree = None

    # Init. libraries
    rna_md = mirmap.if_lib_viennarna.get_model_details(min_loop_size=2, temperature=args.temperature)
    if_spatt = mirmap.if_lib_spatt.Spatt(args.path_libspatt2)

    # Output
    if args.output == "-":
        fout = sys.stdout
    else:
        fout = open(args.output, "wt")
    if args.aggregate:
        if args.output_1to1 == "-":
            fout_1to1 = sys.stdout
        else:
            fout_1to1 = open(args.output_1to1, "wt")
    else:
        fout_1to1 = None

    # Write headers
    if not args.pretty_output:
        fout.write("\t".join(table_columns_annots + table_columns_scores) + "\n")
        if fout_1to1 is not None:
            fout_1to1.write("\t".join(table_columns_annots_1to1 + table_columns_scores) + "\n")

    for (mirna_id, mirna_seq), (transcript_id, transcript_seq) in itertools.product(
        mirnas.items(), transcripts.items()
    ):
        targets = mirmap.target.find_targets_with_seed(
            transcript_seq.upper().replace("U", "T"), mirna_seq.upper().replace("U", "T")
        )

        path_aln = os.path.join(args.path_aln, f"{transcript_id}.fa")
        path_mod = os.path.join(args.path_mod, f"{transcript_id}.mod")
        if not os.path.exists(path_aln):
            path_aln = None
        if not os.path.exists(path_mod):
            path_mod = None

        targets_scores = []
        for itarget, target in enumerate(targets):
            scores = mirmap.scores.calc_scores(
                target,
                rna_md,
                path_aln,
                path_mod,
                tree,
                if_spatt,
                args.path_phylofit,
                args.path_phylop,
            )

            if args.pretty_output:
                out = f"\nTarget #{itarget+1}\n\n" + target.report() + "\n\n" + mirmap.scores.report_scores(scores)
            else:
                out = "\t".join(
                    [
                        mirna_id,
                        transcript_id,
                        str(itarget + 1),
                        str(target.seed_length),
                        str(target.mirna.start),
                        str(target.mirna.end),
                        str(target.seed.start),
                        str(target.seed.end),
                    ]
                    + [format_number(scores[c]) for c in table_columns_scores]
                )
            fout.write(out + "\n")

            targets_scores.append(scores)

        if fout_1to1 is not None and len(targets) > 0:
            seed_lengths = [target.seed_length for target in targets]
            target_agg_scores = mirmap.scores.agg_scores(targets_scores)

            if args.pretty_output:
                out = "\nAggregate scores\n\n" + mirmap.scores.report_scores(target_agg_scores) + "\n"
            else:
                out = "\t".join(
                    [
                        mirna_id,
                        transcript_id,
                        str(len(targets)),
                        str(seed_lengths.count(6)),
                        str(seed_lengths.count(7)),
                    ]
                    + [format_number(target_agg_scores[column]) for column in table_columns_scores]
                )
            fout_1to1.write(out + "\n")


if __name__ == "__main__":
    sys.exit(main())
