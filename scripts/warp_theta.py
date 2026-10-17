# Plots the warp rotation angle theta(r) of the "twist" and "swirl" warp
# types (fractal_type.cc) as a function of the distance r from the window
# center, for any set of parameters:
#
#   twist :  theta(r) = A * (r / r0)^2
#   swirl :  theta(r) = A * exp(-decay * (1 - min(r, r0)/r0))
#
# Both reach A at r = r0 (by default the window corner, r0 = 0.7071), but
# differ in the interior: twist is flat near the center and steep at the
# border, while swirl's exponential falloff (the decay parameter) controls
# how far the twist reaches into the center.
#
# usage:  python scripts/warp_theta.py [-o OUT.svg] [-t TWIST ...] [-s SWIRL ...]
#   -t TWIST    twist curve: <twist_degrees>[,radius]  (default 180,0.7071)
#               (repeatable: each -t adds one curve)
#   -s SWIRL    swirl curve: <twist_degrees>[,radius][,decay]  (default 180,0.7071,3)
#               (repeatable: each -s adds one curve)
#   -o OUT      output file (default: warp-theta.svg)
#
# Without -t/-s the standard set is plotted (twist 180, swirl 180 d=3,
# swirl 180 d=1) -- the comparison used in the README.
#
# Examples:
#   python scripts/warp_theta.py
#   python scripts/warp_theta.py -s 180,0.7071,1 -s 180,0.7071,3 -s 180,0.7071,5
#   python scripts/warp_theta.py -t 90,0.5 -t 180 -s 120,0.5,2 -o compare.svg
#
# Requires matplotlib.

import argparse
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def to_floats(spec):
    try:
        return [float(x) for x in spec.split(",")]
    except ValueError:
        raise argparse.ArgumentTypeError("expected comma-separated numbers, got %r" % spec)

def parse_twist(spec):
    p = to_floats(spec)
    if not 1 <= len(p) <= 2:
        raise argparse.ArgumentTypeError("twist: <twist_degrees>[,radius] (1 or 2 values)")
    a = p[0]
    r0 = p[1] if len(p) > 1 else 0.7071
    if a <= 0 or r0 <= 0:
        raise argparse.ArgumentTypeError("twist: values must be > 0")
    return a, r0

def parse_swirl(spec):
    p = to_floats(spec)
    if not 1 <= len(p) <= 3:
        raise argparse.ArgumentTypeError("swirl: <twist_degrees>[,radius][,decay] (1 to 3 values)")
    a = p[0]
    r0 = p[1] if len(p) > 1 else 0.7071
    decay = p[2] if len(p) > 2 else 3.0
    if a <= 0 or r0 <= 0 or decay <= 0:
        raise argparse.ArgumentTypeError("swirl: values must be > 0")
    return a, r0, decay

def main():
    ap = argparse.ArgumentParser(description="plot the twist/swirl warp rotation angle theta(r) "
                                             "(see the file header for the formulas)")
    ap.add_argument("-t", "--twist", dest="twists", action="append", type=parse_twist,
                    metavar="DEG[,RADIUS]", help="twist curve (repeatable)")
    ap.add_argument("-s", "--swirl", dest="swirls", action="append", type=parse_swirl,
                    metavar="DEG[,RADIUS[,DECAY]]", help="swirl curve (repeatable)")
    ap.add_argument("-o", "--output", default="warp-theta.svg", metavar="OUT.svg")
    args = ap.parse_args()

    twists = args.twists or []
    swirls = args.swirls or []
    if not twists and not swirls:
        twists = [(180.0, 0.7071)]
        swirls = [(180.0, 0.7071, 3.0), (180.0, 0.7071, 1.0)]

    rmax = max([r0 for _, r0 in twists] + [r0 for _, r0, _ in swirls])
    r = np.linspace(0.0, rmax, 1000)

    fig, ax = plt.subplots(figsize=(8, 8))
    curves = [("twist", a, r0, None) for a, r0 in twists]
    curves += [("swirl", a, r0, d) for a, r0, d in swirls]

    for kind, a, r0, decay in curves:
        if kind == "twist":
            theta = a * (r / r0) ** 2
            label = "twist:  theta = A*(r/r0)^2,  A=%g deg, r0=%g" % (a, r0)
        else:
            t = np.minimum(r / r0, 1.0)
            theta = a * np.exp(-decay * (1.0 - t))
            label = "swirl (decay=%g):  theta = A*exp(-d*(1-r/r0)),  A=%g deg, r0=%g" % (decay, a, r0)
        ax.plot(r, theta, lw=2.5, label=label)

    amax = max([a for a, _, _ in swirls] + [a for a, _ in twists])
    ax.set_xlabel(r"distance from window center,  $r$  (window corner = 0.7071)")
    ax.set_ylabel(r"rotation angle  $\theta(r)$  (degrees)")
    ax.set_title("warp rotation angle vs. radius")
    ax.set_xlim(0, rmax)
    ax.set_ylim(0, amax * 1.1 + 1.0)
    ax.grid(True, alpha=0.3)
    if len(curves) > 1:
        ax.legend(loc="upper left", fontsize=9)
    fig.tight_layout()
    fig.savefig(args.output)
    print("saved %s (%d curves)" % (args.output, len(curves)))
    return 0

if __name__ == "__main__":
    sys.exit(main())
