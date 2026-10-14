# Renders a "image plane tour through the ABC cube" animation from the pose
# data dumped by the gen_poses program (gen_poses.cc): the window plane
# swinging through the cube, wireframe cube in perspective view, red/green
# arrows = u/v orientation from the window's top-left corner. Plus, when the
# pose file contains the auxiliary data of gen_poses:
#   CURVE    closed path of the window centers (blue)
#   KEYFRAME the path's keyframe centers (orange rings)
#   HINGE    the hinge edge per frame (orange line, 6sides ease 0/1)
# No path math lives here: the poses are exactly the ones of the C++ flight
# path classes (fractal_type.cc).
#
# usage:  python scripts/cube_tour.py <poses_file> <outdir> [--title TITLE]
#   poses_file  output of "gen_poses -t <6sides|wander> [-p params] -f <frames>"
#   outdir      directory for the PPM frames (default: cube-tour-frames)
#   --title     title text (default: derived from the poses_file name)
#
# Example (the wander tour of the README):
#   gen_poses -t wander -p 10,12345,1.0 -f 256 -o wander-poses.txt
#   python scripts/cube_tour.py wander-poses.txt cube-tour-frames \
#          --title "wander: Catmull-Rom through 10 random keyframes (seed 12345), perspective"
#
# The frames are written to <outdir>/frame-NNN.ppm.
# Encode with e.g.:
#   ffmpeg -framerate 25 -start_number 0 -i <outdir>/frame-%03d.ppm \
#          -c:v libx264 -pix_fmt yuv420p -crf 18 out.mp4

import math
import os
import sys
from PIL import Image, ImageDraw, ImageFont

# ---------------------------------------------------------------------------
# parse the gen_poses output
# ---------------------------------------------------------------------------

def load_poses(path):
    poses = []
    curve = []
    keyframes = []
    hinges = {}
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            t = line.split()
            if t[0] == "POSE":
                v = [float(x) for x in t[2:]]
                poses.append((tuple(v[0:3]), tuple(v[3:6]), tuple(v[6:9])))
            elif t[0] == "CURVE":
                curve.append(tuple(float(x) for x in t[1:4]))
            elif t[0] == "KEYFRAME":
                keyframes.append(tuple(float(x) for x in t[1:4]))
            elif t[0] == "HINGE":
                v = [float(x) for x in t[2:]]
                hinges[int(t[1])] = (tuple(v[0:3]), tuple(v[3:6]))
    return poses, curve, keyframes, hinges

# ---------------------------------------------------------------------------
# window point from a pose (u, v in [0,1] window coordinates)
# ---------------------------------------------------------------------------

CORNS = ((0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0))

def win_pt(corner, u_axis, v_axis, u, v):
    return tuple(corner[i] + u*u_axis[i] + v*v_axis[i] for i in range(3))

# ---------------------------------------------------------------------------
# perspective projection: camera at finite distance, slightly oblique
# ---------------------------------------------------------------------------

def cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])

O = (0.5, 0.5, 0.5)
AZ, EL, R = 40.0, 18.0, 2.6   # azimuth 40 deg (not exactly on a corner), 18 deg elevation
azr, elr = math.radians(AZ), math.radians(EL)
CAM = (O[0] - R * math.cos(elr) * math.cos(azr),
       O[1] - R * math.cos(elr) * math.sin(azr),
       O[2] + R * math.sin(elr))
def norm3(v):
    n = math.sqrt(sum(x * x for x in v))
    return tuple(x / n for x in v)
FWD = norm3(tuple(O[i] - CAM[i] for i in range(3)))
RIGHT = norm3(cross(FWD, (0.0, 0.0, 1.0)))
UP = cross(RIGHT, FWD)

def proj(p):
    v = tuple(p[i] - CAM[i] for i in range(3))
    x = v[0]*RIGHT[0] + v[1]*RIGHT[1] + v[2]*RIGHT[2]
    y = v[0]*UP[0] + v[1]*UP[1] + v[2]*UP[2]
    z = v[0]*FWD[0] + v[1]*FWD[1] + v[2]*FWD[2]
    return (x / z, y / z)

CUBE = [(i, j, k) for i in (0, 1) for j in (0, 1) for k in (0, 1)]
EDGES = []
for i, a in enumerate(CUBE):
    for b in CUBE[i+1:]:
        if sum(x != y for x, y in zip(a, b)) == 1:
            EDGES.append((a, b))

# ---------------------------------------------------------------------------
# scaling / centering (bbox of the 8 cube corners)
# ---------------------------------------------------------------------------

W = H = 512
PC = [proj(c) for c in CUBE]
xs = [p[0] for p in PC]; ys = [p[1] for p in PC]
SC = 0.60 * min(W, H) / max(max(xs) - min(xs), max(ys) - min(ys))
cx0 = (max(xs) + min(xs)) / 2.0
cy0 = (max(ys) + min(ys)) / 2.0
OX = W / 2.0 - cx0 * SC
OY = H / 2.0 + 14 - cy0 * SC

def P2(p):
    x, y = proj(p)
    return (OX + x * SC, OY - y * SC)

def draw_arrow(dr, p, q, color, width=3, head=11):
    (x1, y1), (x2, y2) = p, q
    dr.line((x1, y1, x2, y2), fill=color, width=width)
    ang = math.atan2(y2 - y1, x2 - x1)
    for da in (2.65, -2.65):
        dr.line((x2, y2, x2 + head*math.cos(ang+da), y2 + head*math.sin(ang+da)),
                fill=color, width=width)

def load_font(path, size):
    try:
        return ImageFont.truetype(path, size)
    except OSError:
        return ImageFont.load_default()

font_t = load_font("C:/Windows/Fonts/arial.ttf", 20)
font_s = load_font("C:/Windows/Fonts/arial.ttf", 14)
font_a = load_font("C:/Windows/Fonts/arialbd.ttf", 18)

# ---------------------------------------------------------------------------
# render
# ---------------------------------------------------------------------------

args = [a for a in sys.argv[1:] if a != "--title"]
TITLE_OPT = None
for i, a in enumerate(sys.argv[1:]):
    if a == "--title":
        TITLE_OPT = sys.argv[i + 2]
if len(args) < 2:
    print("usage: python scripts/cube_tour.py <poses_file> <outdir> [--title TITLE]", file=sys.stderr)
    sys.exit(2)
POSES_FILE = args[0]
OUTDIR = args[1]

poses, curve, keyframes, hinges = load_poses(POSES_FILE)
if not poses:
    print("no POSE lines in %s" % POSES_FILE, file=sys.stderr)
    sys.exit(1)
N = len(poses)

TITLE = TITLE_OPT or (os.path.splitext(os.path.basename(POSES_FILE))[0] + ", perspective")
cap = "square = image plane   |   wireframe = cube, perspective view   |   " \
      "red = u, green = v (from top-left corner)"
if curve or keyframes:
    cap += "   |   blue = spline path (center), orange dots = keyframes"
if hinges:
    cap += "   |   orange = hinge edge"

os.makedirs(OUTDIR, exist_ok=True)
for f in range(N):
    pose = poses[f]
    corner, u_axis, v_axis = pose
    img = Image.new('RGB', (W, H), (250, 250, 252))
    dr = ImageDraw.Draw(img)
    tw = dr.textlength(TITLE, font=font_t)
    dr.text(((W - tw) / 2, 14), TITLE, fill=(30, 30, 30), font=font_t)
    cw = dr.textlength(cap, font=font_s)
    dr.text(((W - cw) / 2, H - 26), cap, fill=(110, 110, 110), font=font_s)
    # path of the window centers (closed) and the keyframes, in front of the cube
    if curve:
        pts = [P2(c) for c in curve]
        dr.line(pts + [pts[0]], fill=(40, 120, 220), width=2)
    for c in keyframes:
        x, y = P2(c)
        dr.ellipse((x-5, y-5, x+5, y+5), outline=(230, 120, 0), width=3)
    # cube wireframe: all 12 edges, nothing hidden
    for a, b in EDGES:
        dr.line((P2(a), P2(b)), fill=(70, 70, 78), width=2)
    for c in CUBE:
        x, y = P2(c)
        dr.ellipse((x-3, y-3, x+3, y+3), fill=(70, 70, 78))
    # image plane (quadrilateral), translucent
    quad = [P2(win_pt(corner, u_axis, v_axis, uu, vv)) for (uu, vv) in CORNS]
    ov = Image.new('RGBA', (W, H), (0, 0, 0, 0))
    ovd = ImageDraw.Draw(ov)
    ovd.polygon(quad, fill=(0, 150, 200, 120), outline=(0, 115, 160, 255))
    for q in quad:
        ovd.ellipse((q[0]-4, q[1]-4, q[0]+4, q[1]+4), fill=(0, 115, 160, 255))
    img = Image.alpha_composite(img.convert('RGBA'), ov).convert('RGB')
    dr = ImageDraw.Draw(img)
    # orientation arrows from the window's top-left corner (u,v)=(0,0)
    p00 = P2(corner)
    draw_arrow(dr, p00, P2(win_pt(corner, u_axis, v_axis, 0.42, 0.0)), (215, 60, 60))
    draw_arrow(dr, p00, P2(win_pt(corner, u_axis, v_axis, 0.0, 0.42)), (50, 160, 50))
    # hinge edge in orange (6sides, ease 0/1)
    if f in hinges:
        dr.line((P2(hinges[f][0]), P2(hinges[f][1])), fill=(230, 120, 0), width=4)
    # axis labels: radially outside the corner
    for name, c in (('A', (1, 0, 0)), ('B', (0, 1, 0)), ('C', (0, 0, 1))):
        x, y = P2(c)
        cx, cy = P2(O)
        n = math.hypot(x - cx, y - cy)
        dr.text((x + 20 * (x - cx) / n - 5, y + 20 * (y - cy) / n - 8),
                name, fill=(150, 40, 40), font=font_a)
    img.save(os.path.join(OUTDIR, 'frame-%03d.ppm' % f))
print("%d frames -> %s" % (N, OUTDIR))
