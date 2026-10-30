//! Forced Taylor–Green lattice gas with pressure (see PREREG.md, ADDENDUM.md). Rust 1.81, std only, integer and
//! ternary only: every channel holds a trit 0..2, every diagnostic is an exact integer, json=0.
//! usage: tgwind INIT_TGV.trits GATES.lut
use std::collections::{BTreeMap, BTreeSet};

const K: usize = 27;
const L: usize = 24;
const V: usize = L * L * L;
const T: usize = 480;

#[derive(Clone, Copy)]
struct Vel {
    x: i64,
    y: i64,
    z: i64,
    e: i64,
}

struct Gate {
    q: [usize; 4],
    lut: [u8; 81],
}

fn velocities() -> [Vel; K] {
    let mut v = [Vel { x: 0, y: 0, z: 0, e: 0 }; K];
    let mut i = 0;
    for x in -1..=1i64 {
        for y in -1..=1i64 {
            for z in -1..=1i64 {
                v[i] = Vel { x, y, z, e: x * x + y * y + z * z };
                i += 1;
            }
        }
    }
    v
}

fn decode(mut code: usize) -> [usize; 4] {
    let mut o = [0usize; 4];
    for d in o.iter_mut() {
        *d = code % 3;
        code /= 3;
    }
    o
}

fn encode(d: [u8; 4]) -> usize {
    d[0] as usize + 3 * d[1] as usize + 9 * d[2] as usize + 27 * d[3] as usize
}

/// The repo's gate construction, with ordered maps in place of std::map / std::set.
fn gates(v: &[Vel; K]) -> Vec<Gate> {
    let mut pairs: BTreeMap<(i64, i64, i64, i64), Vec<[usize; 2]>> = BTreeMap::new();
    for a in 0..K {
        for b in a + 1..K {
            pairs.entry((v[a].x + v[b].x, v[a].y + v[b].y, v[a].z + v[b].z, v[a].e + v[b].e)).or_default().push([a, b]);
        }
    }
    let mut quads: BTreeSet<[usize; 4]> = BTreeSet::new();
    for p in pairs.values() {
        for i in 0..p.len() {
            for j in i + 1..p.len() {
                let mut q = [p[i][0], p[i][1], p[j][0], p[j][1]];
                q.sort_unstable();
                if q[0] != q[1] && q[1] != q[2] && q[2] != q[3] {
                    quads.insert(q);
                }
            }
        }
    }
    let mut out = Vec::new();
    for q in quads {
        let mut lut = [0u8; 81];
        let mut cls: BTreeMap<[i64; 5], Vec<usize>> = BTreeMap::new();
        for (mask, l) in lut.iter_mut().enumerate() {
            let d = decode(mask);
            let mut k5 = [0i64; 5];
            for (j, &dj) in d.iter().enumerate() {
                let (w, vel) = (dj as i64, v[q[j]]);
                k5[0] += w;
                k5[1] += w * vel.x;
                k5[2] += w * vel.y;
                k5[3] += w * vel.z;
                k5[4] += w * vel.e;
            }
            cls.entry(k5).or_default().push(mask);
            *l = mask as u8;
        }
        for codes in cls.values() {
            let c: Vec<usize> = codes.iter().copied().filter(|&m| m != 40).collect();
            for pair in c.chunks_exact(2) {
                lut[pair[0]] = pair[1] as u8;
                lut[pair[1]] = pair[0] as u8;
            }
        }
        if lut.iter().enumerate().any(|(m, &l)| l as usize != m) {
            out.push(Gate { q, lut });
        }
    }
    out
}


fn collide(s: &mut [u8], t: usize, g: &[Gate]) {
    for cell in 0..V {
        let base = cell * K;
        for j in 0..32usize {
            let gate = &g[(t * 32 + cell * 17 + j * 37) % g.len()];
            let code = encode([s[base + gate.q[0]], s[base + gate.q[1]], s[base + gate.q[2]], s[base + gate.q[3]]]);
            let d = decode(gate.lut[code] as usize);
            for (k, &qk) in gate.q.iter().enumerate() {
                s[base + qk] = d[k] as u8;
            }
        }
    }
}

fn stream(s: &mut Vec<u8>, v: &[Vel; K]) {
    let mut next = vec![0u8; V * K];
    let w = |a: usize, c: i64| ((a as i64 + c + L as i64) % L as i64) as usize;
    for x in 0..L {
        for y in 0..L {
            for z in 0..L {
                for (ki, vel) in v.iter().enumerate() {
                    let (nx, ny, nz) = (w(x, vel.x), w(y, vel.y), w(z, vel.z));
                    next[((nx * L + ny) * L + nz) * K + ki] = s[((x * L + y) * L + z) * K + ki];
                }
            }
        }
    }
    *s = next;
}

// face channels: index = (x+1)*9 + (y+1)*3 + (z+1)
const PX: usize = 22;
const MX: usize = 4;
const PY: usize = 16;
const MY: usize = 10;

/// One push along an axis: sign +1 moves a particle from the minus face to the plus face.
fn push(s: &mut [u8], cell: usize, plus: usize, minus: usize, sign: i64) -> i64 {
    let (from, to) = if sign > 0 { (minus, plus) } else { (plus, minus) };
    let b = cell * K;
    if s[b + from] >= 1 && s[b + to] <= 1 {
        s[b + from] -= 1;
        s[b + to] += 1;
        1
    } else {
        0
    }
}

#[derive(Clone, Copy, PartialEq)]
enum Arm {
    U,
    F1,
    F2,
    F3,
}

fn force(s: &mut [u8], t: usize, arm: Arm) {
    if arm == Arm::U {
        return;
    }
    for x in 0..L {
        for y in 0..L {
            for z in 0..L {
                let cell = (x * L + y) * L + z;
                if (7 * cell + 13 * t) % 8 != 0 {
                    continue;
                }
                match arm {
                    Arm::F1 => {
                        push(s, cell, PX, MX, sgn_s(x) * sgn_c(y) * sgn_c(z));
                        push(s, cell, PY, MY, -sgn_c(x) * sgn_s(y) * sgn_c(z));
                    }
                    Arm::F2 | Arm::F3 => {
                        let dir = if arm == Arm::F3 && y >= 12 { -1 } else { 1 };
                        let reps = if (t / 20) % 3 == 0 { 2 } else { 1 };
                        for _ in 0..reps {
                            push(s, cell, PX, MX, dir);
                        }
                    }
                    Arm::U => {}
                }
            }
        }
    }
}

struct Row {
    e2: i64,
    n: i64,
    bad: i64,
    px: i64,
    py: i64,
    a: i64,
    wj: i64,
    ck: i64,
    pv: i64,
    b: i64,
    bx: i64,
    by: i64,
    d: i64,
}

fn sgn_s(i: usize) -> i64 {
    if i < 12 {
        1
    } else {
        -1
    }
}

fn sgn_c(i: usize) -> i64 {
    if !(6..18).contains(&i) {
        1
    } else {
        -1
    }
}

fn sgn_2(i: usize) -> i64 {
    if matches!(i % 12, 0 | 1 | 2 | 9 | 10 | 11) {
        1
    } else {
        -1
    }
}

fn measure(s: &[u8], v: &[Vel; K]) -> Row {
    let (mut jx, mut jy, mut jz) = (vec![0i64; V], vec![0i64; V], vec![0i64; V]);
    let mut r = Row { e2: 0, n: 0, bad: 0, px: 0, py: 0, a: 0, wj: 0, ck: 0, pv: 0, b: 0, bx: 0, by: 0, d: 0 };
    let mut p3 = vec![0i64; V];
    let mut rho = vec![0i64; V];
    for cell in 0..V {
        for (ki, vel) in v.iter().enumerate() {
            let n = s[cell * K + ki] as i64;
            if n > 2 {
                r.bad += 1;
            }
            r.n += n;
            rho[cell] += n;
            jx[cell] += n * vel.x;
            jy[cell] += n * vel.y;
            jz[cell] += n * vel.z;
            p3[cell] += n * vel.e;
        }
        r.e2 += p3[cell];
        r.px += jx[cell];
        r.py += jy[cell];
        r.pv += (p3[cell] - 54) * (p3[cell] - 54);
    }
    let id = |x: i64, y: i64, z: i64| {
        let m = |a: i64| a.rem_euclid(L as i64) as usize;
        (m(x) * L + m(y)) * L + m(z)
    };
    let mut blk = vec![[0i64; 3]; 216];
    for x in 0..L as i64 {
        for y in 0..L as i64 {
            for z in 0..L as i64 {
                let wx = (jz[id(x, y + 1, z)] - jz[id(x, y - 1, z)]) - (jy[id(x, y, z + 1)] - jy[id(x, y, z - 1)]);
                let wy = (jx[id(x, y, z + 1)] - jx[id(x, y, z - 1)]) - (jz[id(x + 1, y, z)] - jz[id(x - 1, y, z)]);
                let wz = (jy[id(x + 1, y, z)] - jy[id(x - 1, y, z)]) - (jx[id(x, y + 1, z)] - jx[id(x, y - 1, z)]);
                r.wj += wx * wx + wy * wy + wz * wz;
                let c = id(x, y, z);
                let (ux, uy, uz) = (x as usize, y as usize, z as usize);
                r.a += jx[c] * sgn_s(ux) * sgn_c(uy) * sgn_c(uz) - jy[c] * sgn_c(ux) * sgn_s(uy) * sgn_c(uz);
                r.b += p3[c] * (sgn_2(ux) + sgn_2(uy));
                r.bx += p3[c] * sgn_2(ux);
                r.by += p3[c] * sgn_2(uy);
                r.d += rho[c] * (sgn_2(ux) + sgn_2(uy));
                let bi = (ux / 4 * 6 + uy / 4) * 6 + uz / 4;
                blk[bi][0] += jx[c];
                blk[bi][1] += jy[c];
                blk[bi][2] += jz[c];
            }
        }
    }
    r.ck = blk.iter().map(|m| m[0] * m[0] + m[1] * m[1] + m[2] * m[2]).sum();
    r
}

fn mean(rows: &[Row], f: impl Fn(&Row) -> i64, lo: usize, hi: usize) -> i64 {
    rows[lo..=hi].iter().map(f).sum::<i64>().div_euclid((hi - lo + 1) as i64)
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 3 {
        eprintln!("usage: tgwind INIT_TGV.trits GATES.lut");
        std::process::exit(3);
    }
    let read = |p: &str| std::fs::read(p).unwrap_or_else(|_| std::process::exit(3));
    let (tgv, glut) = (read(&args[1]), read(&args[2]));
    let v = velocities();
    let g = gates(&v);
    let mut mine = Vec::with_capacity(g.len() * 85);
    for gate in &g {
        mine.extend(gate.q.iter().map(|&q| q as u8));
        mine.extend_from_slice(&gate.lut);
    }
    if mine != glut || tgv.len() != V * K || tgv.iter().any(|&n| n > 2) {
        println!("TGW_INPUT|result=INCONCLUSIVE|json=0");
        std::process::exit(3);
    }
    let arms = [(Arm::U, "U"), (Arm::F1, "F1"), (Arm::F2, "F2"), (Arm::F3, "F3")];
    let mut all: Vec<Vec<Row>> = Vec::new();
    for &(arm, name) in &arms {
        let mut s = tgv.clone();
        let mut rows = Vec::with_capacity(T + 1);
        for t in 0..=T {
            let r = measure(&s, &v);
            println!(
                "TGWROW|arm={name}|t={t}|E2={}|N={}|occ_bad={}|PX={}|PY={}|A={}|WJ={}|CK={}|PV={}|B={}|json=0",
                r.e2, r.n, r.bad, r.px, r.py, r.a, r.wj, r.ck, r.pv, r.b
            );
            println!("TGWROW2|arm={name}|t={t}|Bx={}|By={}|D={}|json=0", r.bx, r.by, r.d);
            rows.push(r);
            if t < T {
                collide(&mut s, t, &g);
                force(&mut s, t, arm);
                stream(&mut s, &v);
            }
        }
        all.push(rows);
    }
    let pf = |b: bool| if b { "PASS" } else { "FAIL" };
    let q1 = all.iter().all(|a| a.iter().all(|r| r.bad == 0 && r.e2 == a[0].e2 && r.n == a[0].n));
    let (u, f1, f2, f3) = (&all[0], &all[1], &all[2], &all[3]);
    let a0 = u[0].a;
    let faded = |a: &[Row]| (0..=100).any(|t| a[t].a * 10 <= a0);
    let q2 = faded(u) && u.iter().all(|r| r.px == 0 && r.py == 0);
    let sustained = (100..=T).all(|t| f1[t].a * 10 > a0);
    let q4 = (1..=T).all(|t| f2[t].px >= f2[t - 1].px) && faded(f2);
    let max_step = 2 * 2 * (V as i64 / 8 + 1);
    let ck = |a: &[Row]| mean(a, |r| r.ck, 240, T);
    let q5 = f3.iter().all(|r| r.px.abs() <= max_step) && ck(f3) > 2 * ck(u);
    let bm = |a: &[Row]| mean(a, |r| r.b, 240, T);
    let q7 = bm(f1) > 0 && bm(f1) > 3 * bm(u).abs();
    for (i, &(_, name)) in arms.iter().enumerate() {
        let a = &all[i];
        println!(
            "TGWSUM|arm={name}|A_mean_240_480={}|A_min_100_480={}|CK_mean_240_480={}|WJ_mean_240_480={}|PV_mean_240_480={}|B_mean_240_480={}|PX_480={}|json=0",
            mean(a, |r| r.a, 240, T),
            a[100..=T].iter().map(|r| r.a).min().unwrap_or(0),
            ck(a),
            mean(a, |r| r.wj, 240, T),
            mean(a, |r| r.pv, 240, T),
            bm(a),
            a[T].px
        );
    }
    println!(
        "TGWRESULT|Q1={}|Q2={}|Q3={}|Q4={}|Q5={}|Q7={}|json=0",
        pf(q1),
        pf(q2),
        if sustained { "SUSTAINED" } else { "FADED" },
        pf(q4),
        pf(q5),
        pf(q7)
    );
    let dm = |a: &[Row]| mean(a, |r| r.d, 240, T);
    for (i, &(_, name)) in arms.iter().enumerate() {
        println!("TGWSUM2|arm={name}|D_mean_240_480={}|json=0", dm(&all[i]));
    }
    let q9 = dm(f1) > 0 && dm(f1) > 3 * dm(u).abs();
    println!("TGWRESULT2|Q9={}|json=0", pf(q9));
}
