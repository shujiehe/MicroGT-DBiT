#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Map demultiplexed spatial DNA reads to a bacterial reference with BWA-MEM2.

Input:
  combine.fq from the demultiplexing step.

Main steps:
  bwa-mem2 mem -> samtools view filtering -> samtools sort/index

Output:
  align.sorted.bam
  align.sorted.bam.bai
  mapping.log
"""

import argparse
import logging
import os
import shlex
import subprocess


def file_exists(path):
    return os.path.isfile(path) and os.path.getsize(path) > 0


def run_command(cmd, log=None, check=True, stdout=None, stderr=None):
    if log:
        log.info("RUN: %s", " ".join(cmd))

    proc = subprocess.run(
        cmd,
        stdout=stdout,
        stderr=stderr,
        check=False
    )

    if check and proc.returncode != 0:
        raise RuntimeError(f"Command failed: {' '.join(cmd)}")

    return proc


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Map demultiplexed spatial DNA reads with BWA-MEM2 "
            "and generate a sorted BAM file."
        )
    )

    parser.add_argument(
        "--ref",
        default="bacteria_all.fa",
        help="BWA-MEM2 index prefix. Default: bacteria_all.fa"
    )

    parser.add_argument(
        "--fq",
        required=True,
        help="Input FASTQ file from the demultiplexing step, e.g. combine.fq"
    )

    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory"
    )

    parser.add_argument(
        "--threads",
        type=int,
        default=16,
        help="Number of threads. Default: 16"
    )

    parser.add_argument(
        "--mapq",
        type=int,
        default=10,
        help="Minimum MAPQ for alignment filtering. Default: 10"
    )

    parser.add_argument(
        "--extra_bwa",
        default="",
        help="Additional arguments passed to bwa-mem2 mem. Default: empty"
    )

    parser.add_argument(
        "--rg",
        default="",
        help="Read group string passed to bwa-mem2 mem with -R. '\\t' will be converted to real tabs."
    )

    return parser.parse_args()


def main():
    args = parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    log_fp = os.path.join(args.outdir, "mapping.log")
    logging.basicConfig(
        filename=log_fp,
        filemode="w",
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s"
    )
    log = logging.getLogger("mapping")

    rg = args.rg.replace("\\t", "\t") if args.rg else ""

    log.info("=== Mapping started ===")
    log.info("Reference: %s", args.ref)
    log.info("FASTQ: %s", args.fq)
    log.info("Output directory: %s", args.outdir)
    log.info("Threads: %d", args.threads)
    log.info("MAPQ threshold: %d", args.mapq)

    if args.extra_bwa:
        log.info("Additional bwa-mem2 arguments: %s", args.extra_bwa)

    if rg:
        log.info("Read group: %s", rg)

    bam_unsorted = os.path.join(args.outdir, "align.unsorted.bam")
    bam_filtered = os.path.join(args.outdir, "align.filtered.bam")
    bam_sorted = os.path.join(args.outdir, "align.sorted.bam")

    bwa_cmd = [
        "bwa-mem2",
        "mem",
        "-t",
        str(args.threads),
        "-K",
        "100000000"
    ]

    if args.extra_bwa.strip():
        bwa_cmd += shlex.split(args.extra_bwa.strip())

    if rg:
        bwa_cmd += ["-R", rg]

    bwa_cmd += [args.ref, args.fq]

    samtools_view_cmd = [
        "samtools",
        "view",
        "-@",
        str(args.threads),
        "-b",
        "-F",
        "2308",
        "-q",
        str(args.mapq),
        "-o",
        bam_unsorted
    ]

    log.info("Running bwa-mem2 and initial samtools filtering")
    log.info("RUN: %s | %s", " ".join(bwa_cmd), " ".join(samtools_view_cmd))

    with subprocess.Popen(bwa_cmd, stdout=subprocess.PIPE) as p_bwa, \
            subprocess.Popen(samtools_view_cmd, stdin=p_bwa.stdout) as p_view:

        assert p_bwa.stdout is not None
        p_bwa.stdout.close()

        rc_bwa = p_bwa.wait()
        rc_view = p_view.wait()

        if rc_bwa != 0 or rc_view != 0 or not file_exists(bam_unsorted):
            raise RuntimeError("Mapping or initial filtering failed")

    log.info("Mapping and initial filtering completed")

    os.replace(bam_unsorted, bam_filtered)

    run_command(
        ["samtools", "sort", "-@", str(args.threads), "-o", bam_sorted, bam_filtered],
        log=log
    )

    run_command(
        ["samtools", "index", bam_sorted],
        log=log
    )

    try:
        flagstat = subprocess.check_output(
            ["samtools", "flagstat", bam_sorted],
            text=True
        )
        log.info("\n=== samtools flagstat ===\n%s", flagstat)
    except Exception as exc:
        log.warning("samtools flagstat failed: %s", exc)

    try:
        if os.path.exists(bam_filtered):
            os.remove(bam_filtered)
    except Exception as exc:
        log.warning("Failed to remove intermediate BAM: %s", exc)

    log.info("Output BAM: %s", bam_sorted)
    log.info("Output index: %s.bai", bam_sorted)
    log.info("=== Mapping completed ===")


if __name__ == "__main__":
    main()