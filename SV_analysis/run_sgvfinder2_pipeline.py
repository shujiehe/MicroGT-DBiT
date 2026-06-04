import argparse
from glob import glob
from pathlib import Path

import pandas as pd

from SGVFinder2 import single_file, get_sample_map, work_on_collection
from SGVFinder2.helpers.Bowtie2WrapperSlim import MapPreset


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run SGVFinder2 on single-end R1 FASTQ files."
    )

    parser.add_argument("--db_prefix", required=True, help="SGVFinder2 database prefix.")
    parser.add_argument("--fastq_dir", required=True, help="Directory containing R1 FASTQ files.")
    parser.add_argument("--icra_root", required=True, help="Output directory for per-sample SGVFinder2 files.")
    parser.add_argument("--smp_dir", required=True, help="Output directory for per-sample sample maps.")
    parser.add_argument("--final_dir", required=True, help="Output directory for final SGV results.")

    parser.add_argument("--threads", type=int, default=8, help="Number of threads. Default: 8.")

    parser.add_argument(
        "--sensitivity",
        type=str,
        default="very-sensitive",
        choices=["very-fast", "fast", "sensitive", "very-sensitive"],
        help="Bowtie2 sensitivity preset. Default: very-sensitive."
    )

    parser.add_argument(
        "--report_alignments",
        type=str,
        default="20",
        help='Number of alignments to report, or "all". Default: 20.'
    )

    parser.add_argument("--x_coverage", type=float, default=0.01)
    parser.add_argument("--rate_param", type=int, default=10)
    parser.add_argument("--min_samp_cutoff", type=int, default=2)

    parser.add_argument("--max_spacing", type=int, default=10)
    parser.add_argument("--dels_detect_thresh", type=float, default=0.25)
    parser.add_argument("--real_del_thresh", type=float, default=0.95)
    parser.add_argument("--vsgv_dissim_thresh", type=float, default=0.125)
    parser.add_argument("--dels_cooc_thresh", type=float, default=0.25)
    parser.add_argument("--vsgv_clip_quantile", type=float, default=0.02)
    parser.add_argument("--vsgv_fit_interval", type=float, default=0.95)

    parser.add_argument(
        "--vsgv_fit_method",
        type=str,
        default="betaprime",
        choices=["betaprime", "ncx2"]
    )

    parser.add_argument("--vsgv_dense_perc", type=float, default=85)

    parser.add_argument(
        "--write_csv",
        action="store_true",
        help="Write final SGV tables as CSV files in addition to PKL files."
    )

    return parser.parse_args()


def normalize_report_alignments(value: str):
    if value == "all":
        return value
    return int(value)


def get_sensitivity(value: str):
    mapping = {
        "very-fast": MapPreset.VERY_FAST,
        "fast": MapPreset.FAST,
        "sensitive": MapPreset.SENSITIVE,
        "very-sensitive": MapPreset.VERY_SENSITIVE,
    }
    return mapping[value]


def sample_name_from_r1(r1_path: str) -> str:
    name = Path(r1_path).name

    for suffix in ("_R1.fastq.gz", "_R1.fq.gz"):
        if name.endswith(suffix):
            return name[:-len(suffix)]

    return Path(r1_path).stem


def main():
    args = parse_args()

    Path(args.icra_root).mkdir(parents=True, exist_ok=True)
    Path(args.smp_dir).mkdir(parents=True, exist_ok=True)
    Path(args.final_dir).mkdir(parents=True, exist_ok=True)

    report_alignments = normalize_report_alignments(args.report_alignments)
    senspreset = get_sensitivity(args.sensitivity)

    r1_files = sorted(
        glob(f"{args.fastq_dir}/*_R1.fastq.gz") +
        glob(f"{args.fastq_dir}/*_R1.fq.gz")
    )

    if not r1_files:
        raise RuntimeError(
            f"No R1 FASTQ files found in {args.fastq_dir}. "
            "Expected files ending with _R1.fastq.gz or _R1.fq.gz."
        )

    print(
        f"[INFO] samples={len(r1_files)} "
        f"threads={args.threads} "
        f"sensitivity={args.sensitivity} "
        f"report_alignments={report_alignments}",
        flush=True
    )

    samp_to_map = {}

    for r1 in r1_files:
        sample = sample_name_from_r1(r1)
        print(f"[INFO] processing {sample}", flush=True)

        outdir = Path(args.icra_root) / sample
        outdir.mkdir(parents=True, exist_ok=True)

        _, jsdel_file = single_file(
            fq1=r1,
            fq2=None,
            outfol=str(outdir),
            dbpath=args.db_prefix,
            threads=args.threads,
            senspreset=senspreset,
            report_alns=report_alignments,
        )

        sample_map = get_sample_map(
            delta_fname=jsdel_file,
            lengthdbpath=args.db_prefix + ".dlen",
            x_coverage=args.x_coverage,
            rate_param=args.rate_param,
        )

        smp_file = Path(args.smp_dir) / f"{sample}.smp"
        pd.to_pickle(sample_map, smp_file)
        samp_to_map[sample] = sample_map

    print("[INFO] running collection analysis", flush=True)

    vsgv, dsgv = work_on_collection(
        samp_to_map=samp_to_map,
        max_spacing=args.max_spacing,
        min_samp_cutoff=args.min_samp_cutoff,
        delsdetectthresh=args.dels_detect_thresh,
        real_del_thresh=args.real_del_thresh,
        dels_cooc_thresh=args.dels_cooc_thresh,
        vsgv_dissim_thresh=args.vsgv_dissim_thresh,
        vsgv_clip_quantile=args.vsgv_clip_quantile,
        vsgv_fit_interval=args.vsgv_fit_interval,
        vsgv_fit_method=args.vsgv_fit_method,
        x_coverage=args.x_coverage,
        rate_param=args.rate_param,
        vsgv_dense_perc=args.vsgv_dense_perc,
        browser_path=None,
        taxonomypath=None,
        genepospath=None,
        frames_path=None,
    )

    vsgv_pkl = Path(args.final_dir) / "variable_sgv.pkl"
    dsgv_pkl = Path(args.final_dir) / "deletion_sgv.pkl"

    vsgv.to_pickle(vsgv_pkl)
    dsgv.to_pickle(dsgv_pkl)

    if args.write_csv:
        vsgv.to_csv(Path(args.final_dir) / "variable_sgv.csv")
        dsgv.to_csv(Path(args.final_dir) / "deletion_sgv.csv")

    print(f"[DONE] results written to {args.final_dir}", flush=True)


if __name__ == "__main__":
    main()