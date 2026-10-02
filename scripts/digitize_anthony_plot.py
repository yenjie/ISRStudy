#!/usr/bin/env python3
"""Read C_ISR values off Anthony's tau = 1 - T plot.

This is a screenshot, so the values are approximate: the axis calibration is
taken from the tick labels and the unity line, and each marker is located by
colour.  The reading precision is about +-0.004 in C_ISR; it is stated on every
output and the digitised points are drawn back onto the image for inspection.

Usage:
  scripts/digitize_anthony_plot.py --image <png> --out-csv <csv> --out-check <png>
"""

import argparse
import numpy as np
from PIL import Image, ImageDraw

# matplotlib tab10 colours used in the plot, and which series each one is
SERIES = {
    "KKMCee 5.00.02": (0, 0, 0),
    "Herwig 7.1.5": (31, 119, 180),
    "PYTHIA 8.317": (140, 86, 75),
    "Sherpa 2.2.6": (127, 127, 127),
    "PYTHIA 8.317 Vincia": (23, 190, 207),
}
COLLINEAR = ["Herwig 7.1.5", "PYTHIA 8.317", "Sherpa 2.2.6", "PYTHIA 8.317 Vincia"]


def clusters_1d(idx, gap):
    out = []
    for i in idx:
        if out and i - out[-1][-1] <= gap:
            out[-1].append(i)
        else:
            out.append([i])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--image", required=True)
    ap.add_argument("--out-csv", required=True)
    ap.add_argument("--out-check", required=True)
    ap.add_argument("--nbins", type=int, default=12)
    args = ap.parse_args()

    im = Image.open(args.image).convert("RGB")
    a = np.asarray(im).astype(int)
    H, W, _ = a.shape
    dark = a.sum(axis=2) < 150

    # Frame
    rows = np.where(dark.sum(axis=1) > 0.6 * W)[0]
    cols = np.where(dark.sum(axis=0) > 0.6 * H)[0]
    top, bottom, left, right = rows.min(), rows.max(), cols.min(), cols.max()

    # y calibration: the dashed unity line and the 1.1 major tick on the left axis
    gray = ((abs(a[:, :, 0] - a[:, :, 1]) < 12) & (abs(a[:, :, 1] - a[:, :, 2]) < 12)
            & (a[:, :, 0] > 90) & (a[:, :, 0] < 190))
    g_rows = gray[top + 5:bottom - 5, left + 5:right - 5].sum(axis=1)
    y_unity = top + 5 + int(np.argmax(g_rows))
    tick = []
    for y in range(top + 1, bottom):
        x = left + 1
        while x < left + 40 and dark[y, x]:
            x += 1
        tick.append((y, x - left - 1))
    majors = [y for y, L in tick if L >= 25 and abs(y - y_unity) > 20]
    maj_cl = [np.mean(c) for c in clusters_1d(majors, 3)]
    # nearest major tick above unity is C = 1.1
    y_11 = max(y for y in maj_cl if y < y_unity)
    ppu_y = (y_unity - y_11) / 0.1          # pixels per unit of C

    # x calibration: tick-label centroids below the axis ("0.0" ... "0.4")
    band = dark[bottom + 8:bottom + 60, :]
    colsum = band.sum(axis=0)
    lab_cols = np.where(colsum > 0)[0]
    lab_cl = [c for c in clusters_1d(lab_cols, 14) if 20 < (c[-1] - c[0]) < 90]
    centres = [0.5 * (c[0] + c[-1]) for c in lab_cl]
    centres = [c for c in centres if left - 30 < c < right + 30][:5]
    taus = [0.0, 0.1, 0.2, 0.3, 0.4][:len(centres)]
    slope, x0 = np.polyfit(taus, centres, 1)
    print("frame rows %d-%d cols %d-%d" % (top, bottom, left, right))
    print("y: unity row %d, C=1.1 row %.1f -> %.1f px per unit" % (y_unity, y_11, ppu_y))
    print("x: label centres", [round(c, 1) for c in centres], "-> x = %.1f + %.1f tau" % (x0, slope))

    def to_c(y):
        return 1.0 - (y - y_unity) / ppu_y

    draw = ImageDraw.Draw(im)
    out = ["series,tau,c_isr,reading_precision"]
    summary = {}
    for k in range(args.nbins):
        tau = 0.005 + 0.01 * k
        xc = x0 + slope * tau
        half = int(0.42 * slope * 0.01)           # most of the bin width
        xs = slice(max(left + 2, int(xc - half)), min(right - 2, int(xc + half)))
        ys = slice(top + 3, int(y_unity + 0.13 * ppu_y))   # stop above the legend (C > 0.87)
        sub = a[ys, xs]
        for name, rgb in SERIES.items():
            d = np.abs(sub - np.array(rgb)).sum(axis=2)
            m = d < (45 if name == "KKMCee 5.00.02" else 70 if name != "Sherpa 2.2.6" else 40)
            # marker rows: at least 4 matching pixels wide (error bars are ~2)
            rowsum = m.sum(axis=1)
            mk = np.where(rowsum >= (4 if name == "KKMCee 5.00.02" else 3))[0]
            if len(mk) == 0:
                continue
            cl = clusters_1d(list(mk), 2)
            cl = [c for c in cl if 3 <= len(c) <= 18]
            if not cl:
                continue
            # the marker is the cluster with the largest total pixel count
            best = max(cl, key=lambda c: rowsum[c].sum())
            yrow = ys.start + float(np.average(best, weights=rowsum[best]))
            xcol = xs.start + float(np.average(np.where(m[best].any(axis=0))[0]))
            c = to_c(yrow)
            summary.setdefault(name, []).append((tau, c))
            out.append("%s,%.3f,%.4f,0.004" % (name, tau, c))
            r = 7
            draw.ellipse([xcol - r, yrow - r, xcol + r, yrow + r], outline=(255, 0, 255), width=2)
    # collinear cluster mean per bin
    for k in range(args.nbins):
        tau = 0.005 + 0.01 * k
        vals = [c for n in COLLINEAR for t, c in summary.get(n, []) if abs(t - tau) < 1e-6]
        if vals:
            out.append("collinear models,%.3f,%.4f,%.4f" % (tau, np.mean(vals),
                                                          max(0.004, np.std(vals))))
    with open(args.out_csv, "w") as fh:
        fh.write("\n".join(out) + "\n")
    im.save(args.out_check)
    print("\n%-22s" % "tau" + "".join("%8.3f" % (0.005 + 0.01 * k) for k in range(args.nbins)))
    for name in SERIES:
        row = {round(t, 3): c for t, c in summary.get(name, [])}
        print("%-22s" % name + "".join("%8.3f" % row[round(0.005 + 0.01 * k, 3)]
              if round(0.005 + 0.01 * k, 3) in row else "       -" for k in range(args.nbins)))
    print("wrote", args.out_csv, "and", args.out_check)


if __name__ == "__main__":
    main()
