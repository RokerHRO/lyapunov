# Generates the 6 "image plane tour through the ABC cube" animations shown in
# the README (SixSides section): the window plane swinging around the cube,
# wireframe cube in perspective view, red/green arrows = u/v orientation from
# the window's top-left corner. The poses are exactly the ones of the C++
# fractal type "6sides" (fractal_type.cc).
#
# usage:  python scripts/cube_tour.py [outdir] [--rotate]
#   outdir   directory for the PPM frames (default: cube-tour-frames)
#   --rotate render the variants with in-plane rotation (rotate=1)
#
# The frames are written to <outdir>/ease{0,1,2}/frame-NNN.ppm (256 frames).
# Encode with e.g.:
#   ffmpeg -framerate 25 -start_number 0 -i <outdir>/ease0/frame-%03d.ppm \
#          -c:v libx264 -pix_fmt yuv420p -crf 18 out.mp4

import math
import os
import sys
from PIL import Image, ImageDraw, ImageFont

# ---------------------------------------------------------------------------
# pose math, exact replica of the anonymous namespace in fractal_type.cc
# ---------------------------------------------------------------------------

def sixside_map(seg, th, u, v):
    s, co = math.sin(th), math.cos(th)
    if seg == 0: return (v, 1-(1-u)*co, (1-u)*s)
    if seg == 1: return (1-(1-v)*co, 1-(1-v)*s, 1-u)
    if seg == 2: return (1-u*s, v, 1-u*co)
    if seg == 3: return (1-u, v*co, 1-v*s)
    if seg == 4: return ((1-u)*co, (1-u)*s, 1-v)
    return ((1-v)*s, 1-u, (1-v)*co)

def cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])

def mat2q(ru, rv, n):
    m = [ru[0], ru[1], ru[2], rv[0], rv[1], rv[2], n[0], n[1], n[2]]
    m00, m01, m02, m10, m11, m12, m20, m21, m22 = m
    tr = m00 + m11 + m22
    if tr > 0:
        S = math.sqrt(tr + 1) * 2
        return (0.25*S, (m12-m21)/S, (m20-m02)/S, (m01-m10)/S)
    elif m00 > m11 and m00 > m22:
        S = math.sqrt(1 + m00 - m11 - m22) * 2
        return ((m12-m21)/S, 0.25*S, (m01+m10)/S, (m20+m02)/S)
    elif m11 > m22:
        S = math.sqrt(1 + m11 - m00 - m22) * 2
        return ((m20-m02)/S, (m01+m10)/S, 0.25*S, (m21+m12)/S)
    else:
        S = math.sqrt(1 + m22 - m00 - m11) * 2
        return ((m01-m10)/S, (m20+m02)/S, (m21+m12)/S, 0.25*S)

def qnorm(q):
    n = math.sqrt(sum(x*x for x in q))
    return tuple(x/n for x in q)

def qdot(a, b): return sum(x*y for x, y in zip(a, b))

def qrot(q, v):
    w, x, y, z = q
    vx, vy, vz = v
    tx, ty, tz = 2*(y*vz - z*vy), 2*(z*vx - x*vz), 2*(x*vy - y*vx)
    return (vx + w*tx + (y*tz - z*ty), vy + w*ty + (z*tx - x*tz), vz + w*tz + (x*ty - y*tx))

def h00(t): return (2*t-3)*t*t + 1
def h01(t): return t*t*(3-2*t)
def ht(t): return t*t*(2*t-2)

def cr_quat(qm, q0, q1, qp, t):
    r = tuple(q0[i]*h00(t) + q1[i]*h01(t) + 0.5*(qp[i]-qm[i])*ht(t) for i in range(4))
    return qnorm(r)

def cr_vec(tm, t0, t1, tp, t):
    return tuple(t0[i]*h00(t) + t1[i]*h01(t) + 0.5*(tp[i]-tm[i])*ht(t) for i in range(3))

# keyframes K0..K6 (face poses), as in SixSide::build_cr_path()
QK, TK = [], []
for k in range(7):
    seg = k if k < 6 else 5
    th = math.pi/2 if k == 6 else 0.0
    p0 = sixside_map(seg, th, 0, 0)
    pu = sixside_map(seg, th, 1, 0)
    pv = sixside_map(seg, th, 0, 1)
    ru = tuple(pu[i]-p0[i] for i in range(3))
    rv = tuple(pv[i]-p0[i] for i in range(3))
    TK.append(p0)
    QK.append(qnorm(mat2q(ru, rv, cross(ru, rv))))
for k in range(1, 7):
    if qdot(QK[k-1], QK[k]) < 0: QK[k] = tuple(-x for x in QK[k])
QE = [tuple(-x for x in QK[5]) if qdot(QK[5], QK[0]) < 0 else QK[5],
      tuple(-x for x in QK[1]) if qdot(QK[1], QK[6]) < 0 else QK[1]]
TE = [TK[5], TK[1]]

def pose(t):
    seg = min(int(t), 5)
    tau = t - seg
    qm = QE[0] if seg == 0 else QK[seg-1]
    qp = QE[1] if seg == 5 else QK[seg+2]
    tm = TE[0] if seg == 0 else TK[seg-1]
    tp = TE[1] if seg == 5 else TK[seg+2]
    q = cr_quat(qm, QK[seg], QK[seg+1], qp, tau)
    tr = cr_vec(tm, TK[seg], TK[seg+1], tp, tau)
    return q, tr

# arc-length LUT for the constant-speed reparameterization of the CR path
M = 6 * 512
arclen = [0.0]
prev = pose(0.0)
for i in range(1, M + 1):
    q, tr = pose(6.0 * i / M)
    d = min(1.0, abs(qdot(prev[0], q)))
    arclen.append(arclen[-1] + math.hypot(2.0 * math.acos(d), math.dist(tr, prev[1])))
    prev = (q, tr)
total = arclen[-1]

def pose_cs(t):
    target = t / 6.0 * total
    lo, hi = 0, M
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if arclen[mid] < target: lo = mid
        else: hi = mid
    span = arclen[hi] - arclen[lo]
    frac = (target - arclen[lo]) / span if span > 0 else 0.0
    return pose((lo + frac) * 6.0 / M)

# ---------------------------------------------------------------------------
# window pose, exact replica of SixSide::per_pixel()
# ---------------------------------------------------------------------------

def plane_pt(ease, f, N, u, v, rotate):
    t = 6.0 * (f + 0.5) / N
    if rotate:
        psi = (f + 0.5) * math.pi / N
        ps, pc = math.sin(psi), math.cos(psi)
        du, dv = u - 0.5, v - 0.5
        u = du*pc - dv*ps + 0.5
        v = du*ps + dv*pc + 0.5
    if ease == 2:
        q, tr = pose_cs(t)
        p = qrot(q, (u, v, 0.0))
        return (p[0] + tr[0], p[1] + tr[1], p[2] + tr[2])
    seg = min(int(t), 5)
    tau = t - seg
    theta = (math.pi / 2.0) * (tau * tau * (3.0 - 2.0 * tau) if ease == 1 else tau)
    return sixside_map(seg, theta, u, v)

CORNS = ((0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0))
def corners(ease, f, N, rotate):
    return [plane_pt(ease, f, N, u, v, rotate) for (u, v) in CORNS]

# ---------------------------------------------------------------------------
# perspective projection: camera at finite distance, slightly oblique
# ---------------------------------------------------------------------------

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

# hinge edge per segment (ease 0/1)
HINGE = [((0,1,0),(1,1,0)), ((1,1,0),(1,1,1)), ((1,0,1),(1,1,1)),
         ((0,0,1),(1,0,1)), ((0,0,0),(0,0,1)), ((0,0,0),(0,1,0))]

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

args = [a for a in sys.argv[1:] if a != "--rotate"]
ROTATE = "--rotate" in sys.argv
OUTDIR = args[0] if args else "cube-tour-frames"
TAG = "perspective" + (", rotate" if ROTATE else "")
TITLES = {
    0: "6sides ease=0: constant angular speed (%s)" % TAG,
    1: "6sides ease=1: smoothstep (%s)" % TAG,
    2: "6sides ease=2: quaternion Catmull-Rom (%s)" % TAG,
}
N = 256

for ease in (0, 1, 2):
    outdir = os.path.join(OUTDIR, "ease%d" % ease)
    os.makedirs(outdir, exist_ok=True)
    for f in range(N):
        img = Image.new('RGB', (W, H), (250, 250, 252))
        dr = ImageDraw.Draw(img)
        tw = dr.textlength(TITLES[ease], font=font_t)
        dr.text(((W - tw) / 2, 14), TITLES[ease], fill=(30, 30, 30), font=font_t)
        cap = "square = image plane   |   wireframe = cube, perspective view   |   " \
              "red = u, green = v (from top-left corner)"
        if ease < 2:
            cap += "   |   orange = hinge edge"
        cw = dr.textlength(cap, font=font_s)
        dr.text(((W - cw) / 2, H - 26), cap, fill=(110, 110, 110), font=font_s)
        # cube wireframe: all 12 edges, nothing hidden
        for a, b in EDGES:
            dr.line((P2(a), P2(b)), fill=(70, 70, 78), width=2)
        for c in CUBE:
            x, y = P2(c)
            dr.ellipse((x-3, y-3, x+3, y+3), fill=(70, 70, 78))
        # image plane (quadrilateral), translucent
        quad = [P2(c) for c in corners(ease, f, N, ROTATE)]
        ov = Image.new('RGBA', (W, H), (0, 0, 0, 0))
        ovd = ImageDraw.Draw(ov)
        ovd.polygon(quad, fill=(0, 150, 200, 120), outline=(0, 115, 160, 255))
        for q in quad:
            ovd.ellipse((q[0]-4, q[1]-4, q[0]+4, q[1]+4), fill=(0, 115, 160, 255))
        img = Image.alpha_composite(img.convert('RGBA'), ov).convert('RGB')
        dr = ImageDraw.Draw(img)
        # orientation arrows from the window's top-left corner (u,v)=(0,0)
        p00 = P2(plane_pt(ease, f, N, 0.0, 0.0, ROTATE))
        draw_arrow(dr, p00, P2(plane_pt(ease, f, N, 0.42, 0.0, ROTATE)), (215, 60, 60))
        draw_arrow(dr, p00, P2(plane_pt(ease, f, N, 0.0, 0.42, ROTATE)), (50, 160, 50))
        # hinge edge in orange (ease 0/1)
        seg = min(int(6.0 * (f + 0.5) / N), 5)
        hinge = HINGE[seg] if ease < 2 else None
        if hinge is not None:
            dr.line((P2(hinge[0]), P2(hinge[1])), fill=(230, 120, 0), width=4)
        # axis labels: radially outside the corner
        for name, c in (('A', (1, 0, 0)), ('B', (0, 1, 0)), ('C', (0, 0, 1))):
            x, y = P2(c)
            cx, cy = P2(O)
            n = math.hypot(x - cx, y - cy)
            dr.text((x + 20 * (x - cx) / n - 5, y + 20 * (y - cy) / n - 8),
                    name, fill=(150, 40, 40), font=font_a)
        img.save(os.path.join(outdir, 'frame-%03d.ppm' % f))
    print("ease=%d (rotate=%s): %d frames -> %s" % (ease, ROTATE, N, outdir))
