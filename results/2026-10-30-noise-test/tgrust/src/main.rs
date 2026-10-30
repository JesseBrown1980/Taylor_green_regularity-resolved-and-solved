//! Taylor–Green lattice-gas test (see PREREG.md, ADDENDUM.md). Rust 1.81, std only, integer and
//! ternary only: every channel holds a trit 0..2, every diagnostic is an exact integer, json=0.
//! usage: tgrust INIT_TGV.trits INIT_NULL.trits GATES.lut
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

fn step(s: &mut Vec<u8>, t: usize, g: &[Gate], v: &[Vel; K]) {
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

struct Row {
    e2: i64,
    kj: i64,
    wj: i64,
    a: i64,
    bad: i64,
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

fn measure(s: &[u8], v: &[Vel; K]) -> Row {
    let (mut jx, mut jy, mut jz) = (vec![0i64; V], vec![0i64; V], vec![0i64; V]);
    let mut r = Row { e2: 0, kj: 0, wj: 0, a: 0, bad: 0 };
    for cell in 0..V {
        for (ki, vel) in v.iter().enumerate() {
            let n = s[cell * K + ki] as i64;
            if n > 2 {
                r.bad += 1;
            }
            jx[cell] += n * vel.x;
            jy[cell] += n * vel.y;
            jz[cell] += n * vel.z;
            r.e2 += n * vel.e;
        }
        r.kj += jx[cell] * jx[cell] + jy[cell] * jy[cell] + jz[cell] * jz[cell];
    }
    let id = |x: i64, y: i64, z: i64| {
        let m = |a: i64| a.rem_euclid(L as i64) as usize;
        (m(x) * L + m(y)) * L + m(z)
    };
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
            }
        }
    }
    r
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 4 {
        eprintln!("usage: tgrust INIT_TGV.trits INIT_NULL.trits GATES.lut");
        std::process::exit(3);
    }
    let read = |p: &str| std::fs::read(p).unwrap_or_else(|_| std::process::exit(3));
    let (tgv, nul, glut) = (read(&args[1]), read(&args[2]), read(&args[3]));
    let v = velocities();
    let g = gates(&v);
    let mut mine = Vec::with_capacity(g.len() * 85);
    for gate in &g {
        mine.extend(gate.q.iter().map(|&q| q as u8));
        mine.extend_from_slice(&gate.lut);
    }
    let gates_equal = mine == glut;
    let trits_ok = tgv.len() == V * K && nul.len() == V * K && tgv.iter().chain(nul.iter()).all(|&n| n <= 2);
    if !gates_equal || !trits_ok {
        println!("TGR_INPUT|gates={}|gates_equal_repo={}|trits_ok={}|result=INCONCLUSIVE|json=0", g.len(), gates_equal as i64, trits_ok as i64);
        std::process::exit(3);
    }
    let mut perm: Vec<usize> = (0..V).collect();
    let mut x: u64 = 20261030;
    for i in (1..V).rev() {
        x = (1103515245u64 * x + 12345) % 2147483648;
        let j = (x % (i as u64 + 1)) as usize;
        perm.swap(i, j);
    }
    let mut scr = vec![0u8; V * K];
    for (i, &p) in perm.iter().enumerate() {
        scr[i * K..i * K + K].copy_from_slice(&tgv[p * K..p * K + K]);
    }
    let names = ["TGV", "SCR", "NULL"];
    let mut states = [tgv, scr, nul];
    let mut rows: Vec<Vec<Row>> = Vec::new();
    for s in states.iter_mut() {
        let mut arm = Vec::with_capacity(T + 1);
        for t in 0..=T {
            arm.push(measure(s, &v));
            if t < T {
                step(s, t, &g, &v);
            }
        }
        rows.push(arm);
    }
    for (a, arm) in rows.iter().enumerate() {
        for (t, r) in arm.iter().enumerate() {
            println!("TGROW|arm={}|t={}|E2={}|KJ={}|WJ={}|A={}|occ_bad={}|json=0", names[a], t, r.e2, r.kj, r.wj, r.a, r.bad);
        }
    }
    let p1 = rows.iter().all(|arm| arm.iter().all(|r| r.bad == 0 && r.e2 == arm[0].e2));
    let p2 = rows[2].iter().all(|r| r.kj == 0 && r.wj == 0 && r.a == 0);
    let (a0, as0) = (rows[0][0].a, rows[1][0].a);
    let p3 = a0 > 0 && as0.abs() * 10 <= a0;
    let (mut within, mut above, mut n) = (0i64, 0i64, 0i64);
    for t in 5..=T {
        let (wt, ws) = (rows[0][t].wj, rows[1][t].wj);
        n += 1;
        if (wt - ws).abs() * 10 <= ws {
            within += 1;
        }
        if wt * 10 > ws * 11 {
            above += 1;
        }
    }
    let p4 = if within * 10 >= n * 9 {
        "NOISE"
    } else if above * 10 >= n * 9 {
        "VORTEX"
    } else {
        "MIXED"
    };
    let half = (0..=T).find(|&t| rows[0][t].a * 2 <= a0).map_or(-1, |t| t as i64);
    let pf = |b: bool| if b { "PASS" } else { "FAIL" };
    println!("TGSUM|ticks_5_480={n}|within_10pct={within}|tgv_above_110pct={above}|P4={p4}|json=0");
    println!(
        "TGSUM|A_TGV0={a0}|A_SCR0={as0}|A_TGV_half_tick={half}|WJ_TGV0={}|WJ_TGV1={}|WJ_SCR0={}|WJ_SCR1={}|json=0",
        rows[0][0].wj, rows[0][1].wj, rows[1][0].wj, rows[1][1].wj
    );
    println!("TGRESULT|P1={}|P2={}|P3={}|P4={p4}|json=0", pf(p1), pf(p2), pf(p3));
}
