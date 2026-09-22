#!/usr/bin/env python3
"""Fill the KKMC measurements into the ISR model comparison slide source.

The deck declares one \\newcommand per KKMC number with a placeholder body.  This
script replaces each placeholder with the measured value taken from the macro
outputs, so the slides can never drift from the files they were made from.

Usage:
  scripts/fill_kkmc_slide_numbers.py \\
      --results /data2/yjlee/ISRsample/kkmc_1M_20260921/results \\
      --runlog  /raid5/data/yjlee/ISR/KKMCee/ffbench/ISRrun_ISR_ON/kkmc_run.log \\
      --slides  /raid5/data/yjlee/ISR/overleaf/slides/20260921-isr_model_comparison.tex
"""

import argparse
import csv
import os
import re
import sys


def read_summary(path):
    """Return {(sample, state): row} from isr_model_summary.csv."""
    out = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            out[(row["sample"], row["isr_state"])] = row
    return out


def read_endpoint(path):
    out = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            out[row["sample"]] = row
    return out


def read_ifi(path):
    out = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            out[row["quantity"]] = (float(row["value"]), float(row["error"]))
    return out


def read_xsec(path):
    """Pull 'xSecPb = <value> +- <error>' out of the KKMC run log, in nb."""
    pat = re.compile(r"xSecPb\s*=\s*([0-9.Ee+-]+)\s*\+-\s*([0-9.Ee+-]+)")
    with open(path, errors="ignore") as fh:
        for line in fh:
            m = pat.search(line)
            if m:
                return float(m.group(1)) / 1000.0, float(m.group(2)) / 1000.0
    raise SystemExit("could not find xSecPb in %s" % path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True)
    ap.add_argument("--runlog", required=True)
    ap.add_argument("--slides", required=True)
    ap.add_argument("--label", default="KKMC 4.30 CEEX")
    args = ap.parse_args()

    summary = read_summary(os.path.join(args.results, "isr_model_summary.csv"))
    endpoint = read_endpoint(os.path.join(args.results, "isr_model_endpoint.csv"))
    ifi = read_ifi(os.path.join(args.results, "isr_model_kkmc_ifi.csv"))
    xsec, xsec_err = read_xsec(args.runlog)

    on = summary[(args.label, "ON")]
    off = summary[(args.label, "OFF")]
    ep = endpoint[args.label]

    mvis_on = float(on["mean_m_vis"])
    mvis_off = float(off["mean_m_vis"])
    ifi_int, ifi_int_err = ifi["integrated_ratio"]
    ifi_last, ifi_last_err = ifi["lastbin_ratio"]

    values = {
        "kkEisrOn": "%.3f" % float(on["mean_e_isr_generator_truth"]),
        "kkEbcOn": "%.3f" % float(on["mean_e_beam_collinear_gamma"]),
        "kkEbcOff": "%.4f" % float(off["mean_e_beam_collinear_gamma"]),
        "kkMvisOff": "%.3f" % mvis_off,
        "kkMvisOn": "%.3f" % mvis_on,
        "kkDeltaMvis": "%+.3f" % (mvis_on - mvis_off),
        "kkTOff": "%.5f" % float(ep["mean_T_off"]),
        "kkTOn": "%.5f" % float(ep["mean_T_on"]),
        "kkDeltaT": "%+.5f" % float(ep["delta_mean_T"]),
        "kkLastBin": "%.3f" % float(ep["lastbin_cisr"]),
        "kkLastBinErr": "%.3f" % float(ep["lastbin_err"]),
        # These two are used outside math mode in the deck, so they carry
        # their own delimiters.
        "kkIfiMean": "$%.4f \\pm %.4f$" % (ifi_int, ifi_int_err),
        "kkIfiLast": "$%.3f \\pm %.3f$" % (ifi_last, ifi_last_err),
        "kkXsec": "%.2f" % xsec,
        "kkXsecErr": "%.2f" % xsec_err,
    }

    with open(args.slides) as fh:
        tex = fh.read()

    missing = []
    for name, value in values.items():
        pat = re.compile(r"(\\newcommand\{\\%s\}\{)[^}]*(\})" % name)
        tex, n = pat.subn(lambda m: m.group(1) + value + m.group(2), tex)
        if n == 0:
            missing.append(name)

    if missing:
        print("warning: no placeholder found for: %s" % ", ".join(missing),
              file=sys.stderr)

    with open(args.slides, "w") as fh:
        fh.write(tex)

    for name in sorted(values):
        print("  \\%-14s = %s" % (name, values[name]))
    print("[done] filled %s" % args.slides)


if __name__ == "__main__":
    main()
