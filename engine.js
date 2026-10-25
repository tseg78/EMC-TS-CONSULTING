/* EMC TS CONSULTING — Déclinaison BCI EMCTS V1 — moteur de calcul
 * Susceptibilité conduite RF (BCI) — modèle de ligne multiconductrice (MTL).
 *
 * Principes :
 *  - Faisceau et câblage interne = lignes multiconductrices sans pertes en milieu homogène
 *    (air, εr effectif réglable). Matrices L′ et C′ calculées depuis la géométrie
 *    (méthode des images, fils fins), avec inductances et capacités mutuelles.
 *    En milieu homogène L′·C′ = μ0·ε·I : la solution par paramètres de chaîne est exacte.
 *  - Pince d'injection = source de tension SÉRIE identique sur chaque conducteur à la
 *    position de la pince (transformateur entourant tout le faisceau).
 *  - Résolution : analyse nodale modifiée (MNA) complexe ; les liaisons idéales (0 Ω)
 *    sont traitées par fusion de nœuds, sans impédance numérique fictive.
 *  - Aucune valeur n'est imposée en dehors de la configuration fournie.
 */
(function (global) {
  'use strict';

  const C0 = 299792458;
  // Nombre maximal de conducteurs : le mode détaillé (matrices pleines) devient lent au-delà de 40.
  const MAX_WIRES = { approx: 60, detailed: 40 };
  const MU0 = 4e-7 * Math.PI;
  const EPS0 = 1 / (MU0 * C0 * C0);

  /* ---------- Nombres complexes ---------- */
  const cx = (r = 0, i = 0) => ({ r, i });
  const add = (a, b) => cx(a.r + b.r, a.i + b.i);
  const sub = (a, b) => cx(a.r - b.r, a.i - b.i);
  const mul = (a, b) => cx(a.r * b.r - a.i * b.i, a.r * b.i + a.i * b.r);
  const div = (a, b) => { const d = b.r * b.r + b.i * b.i; if (!(d > 0)) return null; return cx((a.r * b.r + a.i * b.i) / d, (a.i * b.r - a.r * b.i) / d); };
  const abs = a => Math.hypot(a.r, a.i);
  const scl = (a, s) => cx(a.r * s, a.i * s);
  const inv = a => div(cx(1, 0), a);
  const csin = z => cx(Math.sin(z.r) * Math.cosh(z.i), Math.cos(z.r) * Math.sinh(z.i));
  const ccos = z => cx(Math.cos(z.r) * Math.cosh(z.i), -Math.sin(z.r) * Math.sinh(z.i));
  const cexp = z => { const m = Math.exp(z.r); return cx(m * Math.cos(z.i), m * Math.sin(z.i)); };
  // exp(z) − 1 sans perte de précision pour |z| petit
  const cexpm1 = z => {
    if (abs(z) > 1e-3) return sub(cexp(z), cx(1, 0));
    let t = z, s = z; for (let k = 2; k <= 6; k++) { t = scl(mul(t, z), 1 / k); s = add(s, t); } return s;
  };
  // racine carrée complexe, partie réelle positive
  const csqrt = z => { const m = abs(z); let r = Math.sqrt((m + z.r) / 2), i = Math.sqrt(Math.max(0, (m - z.r) / 2)); if (z.i < 0) i = -i; return cx(r, i); };

  /* ---------- Résistance linéique d'un conducteur en cuivre (effet de peau) ---------- */
  const SIGMA_CU = 5.8e7; // S/m, cuivre recuit
  function rDC(rw) { return rw > 0 ? 1 / (SIGMA_CU * Math.PI * rw * rw) : 0; }
  function skinDepth(f) { return 1 / Math.sqrt(Math.PI * f * MU0 * SIGMA_CU); }
  // R′(f) = R′dc·√(1 + (r/2δ)²) : R′dc en basse fréquence, 1/(2π r δ σ) en haute fréquence.
  function rPrime(f, rw) { if (!(rw > 0)) return 0; const x = rw / (2 * skinDepth(f)); return rDC(rw) * Math.sqrt(1 + x * x); }

  /* ---------- Unités ---------- */
  const unitL = u => (u === 'mH' ? 1e-3 : u === 'uH' ? 1e-6 : 1e-9);
  const unitC = u => (u === 'uF' ? 1e-6 : u === 'nF' ? 1e-9 : 1e-12);
  const num = v => (v === '' || v === null || v === undefined ? NaN : Number(v));
  const pos = v => Number.isFinite(num(v)) && num(v) > 0;
  const nonneg = v => Number.isFinite(num(v)) && num(v) >= 0;

  /* ---------- Algèbre réelle ---------- */
  function invertReal(M) {
    const n = M.length, A = M.map((row, i) => [...row, ...Array.from({ length: n }, (_, j) => (i === j ? 1 : 0))]);
    for (let k = 0; k < n; k++) {
      let p = k;
      for (let i = k + 1; i < n; i++) if (Math.abs(A[i][k]) > Math.abs(A[p][k])) p = i;
      if (Math.abs(A[p][k]) < 1e-300) return null;
      [A[k], A[p]] = [A[p], A[k]];
      const d = A[k][k];
      for (let j = 0; j < 2 * n; j++) A[k][j] /= d;
      for (let i = 0; i < n; i++) if (i !== k) { const q = A[i][k]; if (q) for (let j = 0; j < 2 * n; j++) A[i][j] -= q * A[k][j]; }
    }
    return A.map(row => row.slice(n));
  }

  // Valeurs et vecteurs propres d'une matrice réelle symétrique (méthode de Jacobi). T : vecteurs propres en colonnes.
  function jacobiEig(A0) {
    const n = A0.length, A = A0.map(r => r.slice()), T = Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => (i === j ? 1 : 0)));
    for (let sweep = 0; sweep < 100; sweep++) {
      let off = 0; for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) off += A[i][j] * A[i][j];
      if (off < 1e-30 * (1 + A.reduce((s, r, i) => s + r[i] * r[i], 0))) break;
      for (let p = 0; p < n; p++) for (let q = p + 1; q < n; q++) {
        if (Math.abs(A[p][q]) < 1e-300) continue;
        const th = (A[q][q] - A[p][p]) / (2 * A[p][q]);
        const t = Math.sign(th || 1) / (Math.abs(th) + Math.sqrt(th * th + 1)), c = 1 / Math.sqrt(t * t + 1), sn = t * c;
        for (let k = 0; k < n; k++) { const akp = A[k][p], akq = A[k][q]; A[k][p] = c * akp - sn * akq; A[k][q] = sn * akp + c * akq; }
        for (let k = 0; k < n; k++) { const apk = A[p][k], aqk = A[q][k]; A[p][k] = c * apk - sn * aqk; A[q][k] = sn * apk + c * aqk; }
        for (let k = 0; k < n; k++) { const tkp = T[k][p], tkq = T[k][q]; T[k][p] = c * tkp - sn * tkq; T[k][q] = sn * tkp + c * tkq; }
      }
    }
    return { vals: A.map((r, i) => r[i]), T };
  }

  /* ---------- Géométrie ---------- */
  // Conducteurs répartis sur un anneau dans le faisceau ; centre du faisceau à la hauteur h.
  function ringPositions(n, bundleD, isoD, h) {
    const warnings = [];
    if (n === 1) return { pts: [{ x: 0, y: h }], ring: 0, warnings };
    const aMin = isoD / (2 * Math.sin(Math.PI / n));
    let a = (bundleD - isoD) / 2;
    if (a < aMin) {
      warnings.push(`Diamètre insuffisant pour ${n} conducteurs isolés de ${(isoD * 1e3).toFixed(2)} mm sur un anneau : rayon d'anneau porté à ${(aMin * 1e3).toFixed(2)} mm (faisceau équivalent ${((2 * aMin + isoD) * 1e3).toFixed(1)} mm).`);
      a = aMin;
    }
    const pts = Array.from({ length: n }, (_, k) => ({ x: a * Math.cos(2 * Math.PI * k / n), y: h + a * Math.sin(2 * Math.PI * k / n) }));
    return { pts, ring: a, warnings };
  }
  // Diamètre minimal du cercle contenant N cercles égaux de diamètre 1 (empilements optimaux connus).
  const PACKING = [0, 1, 2, 2.1547, 2.4142, 2.7013, 3, 3, 3.3048, 3.6131, 3.8130, 3.9238, 4.0296, 4.2361, 4.3284, 4.5214, 4.6154];
  function packingRatio(n) { return n < PACKING.length ? PACKING[n] : 1.1 * Math.sqrt(n / 0.9069); }
  // Diamètre extérieur d'un toron de n conducteurs isolés de diamètre isoD (sans gaine).
  function autoBundleD(n, isoD) { return packingRatio(n) * isoD; }

  // Répartition en couches concentriques dans un toron de diamètre bundleD.
  function bundlePositions(n, bundleD, isoD, h) {
    const warnings = [], pts = [];
    if (n === 1) return { pts: [{ x: 0, y: h }], ring: 0, warnings };
    let left = n, rho = Math.max((bundleD - isoD) / 2, 0), layer = 0;
    while (left > 0) {
      if (rho < isoD / 2 || left === 1) {
        if (left === 1) { pts.push({ x: 0, y: h }); left = 0; break; }
        rho = Math.max(rho, isoD / (2 * Math.sin(Math.PI / left)));
      }
      const cap = Math.max(1, Math.floor(Math.PI / Math.asin(Math.min(1, isoD / (2 * rho))) + 0.01));
      let m = left <= cap ? left : cap;
      if (left > cap && rho - isoD < isoD / 2 && !(left === cap + 1 && rho >= isoD * 0.999)) {
        // Dernière couche possible : tous les conducteurs restants sur cet anneau, élargi si nécessaire.
        m = left; rho = Math.max(rho, isoD / (2 * Math.sin(Math.PI / left)));
      }
      // Décalage angulaire de la couche choisi pour maximiser l'écart avec les couches déjà placées.
      const ringAt = off => Array.from({ length: m }, (_, k) => ({ x: rho * Math.cos(off + 2 * Math.PI * k / m), y: h + rho * Math.sin(off + 2 * Math.PI * k / m) }));
      let best = ringAt(0), bestD = -1;
      if (layer > 0) for (let t = 0; t < 36; t++) {
        const cand = ringAt(t / 36 * 2 * Math.PI / m);
        const d = Math.min(...cand.flatMap(p => pts.map(q => Math.hypot(p.x - q.x, p.y - q.y))));
        if (d > bestD) { bestD = d; best = cand; }
      }
      pts.push(...best);
      left -= m; rho -= isoD; layer++;
    }
    const minD = autoBundleD(n, isoD);
    if (bundleD < minD * 0.999) warnings.push(`Le diamètre saisi (${(bundleD * 1e3).toFixed(2)} mm) est inférieur au diamètre minimal d'un toron de ${n} conducteurs isolés de ${(isoD * 1e3).toFixed(2)} mm (${(minD * 1e3).toFixed(2)} mm).`);
    return { pts, ring: (bundleD - isoD) / 2, warnings };
  }
  function rowPositions(n, pitch, h) {
    return { pts: Array.from({ length: n }, (_, k) => ({ x: (k - (n - 1) / 2) * pitch, y: h })), ring: 0, warnings: [] };
  }

  // Matrice d'inductance linéique (fils fins au-dessus d'un plan parfait, méthode des images).
  function inductanceMatrix(pts, rw) {
    const n = pts.length, L = Array.from({ length: n }, () => Array(n).fill(0));
    for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) {
      if (i === j) L[i][i] = MU0 / (2 * Math.PI) * Math.log(2 * pts[i].y / rw);
      else {
        const dx = pts[i].x - pts[j].x, dy = pts[i].y - pts[j].y, d2 = dx * dx + dy * dy;
        L[i][j] = MU0 / (4 * Math.PI) * Math.log(1 + 4 * pts[i].y * pts[j].y / d2);
      }
    }
    return L;
  }

  // Paramètres complets d'une ligne homogène.
  function lineParams(geo, rw, epsr) {
    const errors = [];
    for (const p of geo.pts) if (p.y <= rw) errors.push('Un conducteur touche ou traverse le plan de référence : hauteur insuffisante.');
    if (errors.length) return { errors, warnings: geo.warnings };
    const L = inductanceMatrix(geo.pts, rw);
    const Linv = invertReal(L);
    if (!Linv) return { errors: ['Matrice d\'inductance singulière.'], warnings: geo.warnings };
    const v = C0 / Math.sqrt(epsr);
    const C = Linv.map(r => r.map(x => x * MU0 * EPS0 * epsr));
    const Yc = Linv.map(r => r.map(x => x / v));
    const sumC = C.flat().reduce((a, b) => a + b, 0);
    const sumYc = Yc.flat().reduce((a, b) => a + b, 0);
    return { L, C, Yc, v, epsr, n: L.length, Ccm: sumC, Zcm: 1 / sumYc, eig: jacobiEig(C), errors, warnings: geo.warnings, pts: geo.pts, rw };
  }

  /* ---------- Circuit MNA ---------- */
  class Circuit {
    constructor() { this.count = 1; this.parent = [0]; this.elems = []; this.srcs = []; }
    node() { this.parent.push(this.count); return this.count++; }
    nodes(n) { return Array.from({ length: n }, () => this.node()); }
    find(a) { while (this.parent[a] !== a) { this.parent[a] = this.parent[this.parent[a]]; a = this.parent[a]; } return a; }
    short(a, b) {
      a = this.find(a); b = this.find(b); if (a === b) return;
      if (a === 0) this.parent[b] = 0; else if (b === 0) this.parent[a] = 0; else this.parent[b] = a;
    }
    // Impédance entre deux nœuds ; null = circuit ouvert ; |Z| = 0 = liaison idéale.
    z(a, b, Z) {
      if (!Z) return;
      if (Z.r === 0 && Z.i === 0) { this.short(a, b); return; }
      this.elems.push({ t: 'y', a, b, Y: inv(Z) });
    }
    // Ligne avec pertes : résistance linéique Rp identique sur chaque conducteur (plan de référence parfait).
    // Z′ = Rp·I + jωL′, Y′ = jωC′, L′ = C′⁻¹/v² : tout est diagonal dans la base propre de C′ (C′ = T Λ Tᵀ).
    // Mode m : z = Rp + jω/(v²λm), y = jωλm, γ = √(zy), Y0 = y/γ ; Y11 = T diag(Y0·coth γl) Tᵀ, Y12 = −T diag(Y0·csch γl) Tᵀ.
    line(p1, p2, length, P, w, Rp = 0) {
      if (Rp > 0) {
        const { vals, T } = P.eig, n = vals.length, a = [], b = [];
        for (let m = 0; m < n; m++) {
          const z = cx(Rp, w / (P.v * P.v * vals[m])), y = cx(0, w * vals[m]);
          const g = csqrt(mul(z, y)), Y0 = div(y, g), gl = scl(g, length);
          const den = scl(cexpm1(scl(gl, -2)), -1); // 1 − e^(−2γl)
          const coth = div(add(cx(1, 0), cexp(scl(gl, -2))), den), csch = div(scl(cexp(scl(gl, -1)), 2), den);
          a.push(mul(Y0, coth)); b.push(scl(mul(Y0, csch), -1));
        }
        const A = [], B = [];
        for (let i = 0; i < n; i++) { A.push([]); B.push([]); for (let j = 0; j < n; j++) {
          let sa = cx(0, 0), sb = cx(0, 0);
          for (let m = 0; m < n; m++) { const t = T[i][m] * T[j][m]; if (t) { sa = add(sa, scl(a[m], t)); sb = add(sb, scl(b[m], t)); } }
          A[i].push(sa); B[i].push(sb);
        } }
        this.elems.push({ t: 'mtlc', p1, p2, A, B });
        return;
      }
      const theta = cx(w * length / P.v, 0);
      let s = csin(theta);
      if (abs(s) < 1e-12) s = csin(cx(theta.r, -1e-9)); // résonance exacte : régularisation infinitésimale
      const c = ccos(theta);
      this.elems.push({ t: 'mtl', p1, p2, a: mul(cx(0, -1), div(c, s)), b: mul(cx(0, 1), inv(s)), Yc: P.Yc });
    }
    vsrc(a, b) { this.srcs.push({ a, b }); return this.srcs.length - 1; }

    // Sources : V(b) − V(a) = E − Zp · Σ I_k  (même f.é.m. sur chaque conducteur, impédance commune Zp).
    solve(E, Zp) {
      const map = new Map(); let n = 0;
      for (let a = 1; a < this.count; a++) { const r = this.find(a); if (r !== 0 && !map.has(r)) map.set(r, n++); }
      const idx = a => { const r = this.find(a); return r === 0 ? -1 : map.get(r); };
      const ns = this.srcs.length, dim = n + ns;
      const Mr = new Float64Array(dim * dim), Mi = new Float64Array(dim * dim), Br = new Float64Array(dim), Bi = new Float64Array(dim);
      const st = (i, j, y) => { if (i < 0 || j < 0) return; Mr[i * dim + j] += y.r; Mi[i * dim + j] += y.i; };
      for (const e of this.elems) {
        if (e.t === 'y') {
          const i = idx(e.a), j = idx(e.b), y = e.Y;
          if (i === j) continue;
          st(i, i, y); st(j, j, y); st(i, j, scl(y, -1)); st(j, i, scl(y, -1));
        } else if (e.t === 'mtlc') {
          const N = e.p1.length, i1 = e.p1.map(idx), i2 = e.p2.map(idx);
          for (let i = 0; i < N; i++) for (let j = 0; j < N; j++) {
            st(i1[i], i1[j], e.A[i][j]); st(i2[i], i2[j], e.A[i][j]); st(i1[i], i2[j], e.B[i][j]); st(i2[i], i1[j], e.B[i][j]);
          }
        } else {
          const N = e.p1.length;
          const i1 = e.p1.map(idx), i2 = e.p2.map(idx);
          for (let i = 0; i < N; i++) for (let j = 0; j < N; j++) {
            const y = e.Yc[i][j], A = scl(e.a, y), B = scl(e.b, y);
            st(i1[i], i1[j], A); st(i2[i], i2[j], A); st(i1[i], i2[j], B); st(i2[i], i1[j], B);
          }
        }
      }
      this.srcs.forEach((s, k) => {
        const r = n + k, a = idx(s.a), b = idx(s.b);
        if (a >= 0) { Mr[a * dim + r] += 1; Mr[r * dim + a] -= 1; }
        if (b >= 0) { Mr[b * dim + r] -= 1; Mr[r * dim + b] += 1; }
        for (let j = 0; j < ns; j++) { Mr[r * dim + n + j] += Zp.r; Mi[r * dim + n + j] += Zp.i; }
        Br[r] = E.r; Bi[r] = E.i;
      });
      // Élimination de Gauss complexe avec pivot partiel.
      for (let k = 0; k < dim; k++) {
        let p = k, best = Mr[k * dim + k] ** 2 + Mi[k * dim + k] ** 2;
        for (let i = k + 1; i < dim; i++) { const m = Mr[i * dim + k] ** 2 + Mi[i * dim + k] ** 2; if (m > best) { best = m; p = i; } }
        if (!(best > 1e-300)) return null;
        if (p !== k) {
          for (let j = 0; j < dim; j++) {
            let t = Mr[k * dim + j]; Mr[k * dim + j] = Mr[p * dim + j]; Mr[p * dim + j] = t;
            t = Mi[k * dim + j]; Mi[k * dim + j] = Mi[p * dim + j]; Mi[p * dim + j] = t;
          }
          let t = Br[k]; Br[k] = Br[p]; Br[p] = t; t = Bi[k]; Bi[k] = Bi[p]; Bi[p] = t;
        }
        const pr = Mr[k * dim + k], pi = Mi[k * dim + k];
        for (let i = k + 1; i < dim; i++) {
          const ar = Mr[i * dim + k], ai = Mi[i * dim + k];
          if (ar === 0 && ai === 0) continue;
          const qr = (ar * pr + ai * pi) / best, qi = (ai * pr - ar * pi) / best;
          for (let j = k; j < dim; j++) {
            const xr = Mr[k * dim + j], xi = Mi[k * dim + j];
            Mr[i * dim + j] -= qr * xr - qi * xi; Mi[i * dim + j] -= qr * xi + qi * xr;
          }
          Br[i] -= qr * Br[k] - qi * Bi[k]; Bi[i] -= qr * Bi[k] + qi * Br[k];
        }
      }
      const Xr = new Float64Array(dim), Xi = new Float64Array(dim);
      for (let i = dim - 1; i >= 0; i--) {
        let sr = Br[i], si = Bi[i];
        for (let j = i + 1; j < dim; j++) { const ar = Mr[i * dim + j], ai = Mi[i * dim + j]; sr -= ar * Xr[j] - ai * Xi[j]; si -= ar * Xi[j] + ai * Xr[j]; }
        const dr = Mr[i * dim + i], di = Mi[i * dim + i], dd = dr * dr + di * di;
        Xr[i] = (sr * dr + si * di) / dd; Xi[i] = (si * dr - sr * di) / dd;
      }
      return {
        V: a => { const i = idx(a); return i < 0 ? cx(0, 0) : cx(Xr[i], Xi[i]); },
        I: k => cx(Xr[n + k], Xi[n + k]),
        dim
      };
    }
  }

  /* ---------- Composants ---------- */
  function capZ(w, C, esr, esl) { if (!(C > 0)) return null; return cx(esr || 0, w * (esl || 0) - 1 / (w * C)); }
  function indZ(w, L, dcr, cpar) {
    if (!(L > 0)) return null;
    const zs = cx(dcr || 0, w * L);
    if (cpar > 0) return inv(add(inv(zs), cx(0, w * cpar)));
    return zs;
  }
  // Terminaison côté charge d'un signal / d'une alimentation.
  // Résistance réelle : R avec inductance série Ls (nH) et capacité parallèle Cp (pF) facultatives (case vide = non modélisée).
  function realR(R, lsNh, cpPf, w) {
    let z = cx(R, w * (num(lsNh) || 0) * 1e-9);
    const cp = (num(cpPf) || 0) * 1e-12;
    if (cp > 0) z = inv(add(inv(z), cx(0, w * cp)));
    return z;
  }
  function termZ(e, w) {
    const R = num(e.r), C = num(e.c) * unitC(e.cunit), L = num(e.l) * unitL(e.lunit);
    switch (e.type) {
      case 'R': return realR(R, e.ls, e.cp, w);
      case 'Z': return cx(num(e.z), 0);
      case 'C': return capZ(w, C);
      case 'RCs': return cx(R, -1 / (w * C));
      case 'RCp': return inv(cx(1 / R, w * C));
      case 'LC': return cx(0, w * L - 1 / (w * C));
      default: return null; // ouvert
    }
  }
  // Charge côté carte.
  function loadZ(d, w) {
    const R = num(d.loadR), L = num(d.loadL) * unitL(d.loadLUnit), C = num(d.loadC) * unitC(d.loadCUnit);
    switch (d.loadType) {
      case 'R': return realR(R, d.loadLs, d.loadCp, w);
      case 'Z': return cx(R, 0);
      case 'C': return capZ(w, C);
      case 'L': return cx(0, w * L);
      case 'RCs': return cx(R, -1 / (w * C));
      case 'RCp': return inv(cx(1 / R, w * C));
      case 'LC': return cx(0, w * L - 1 / (w * C));
      default: return null;
    }
  }
  // Filtre série/shunt entre inN et outN, branche shunt vers refNode.
  function stampFilter(ck, inN, outN, d, refNode, w) {
    const R = num(d.r) || 0, L = (num(d.l) || 0) * unitL(d.lunit), C = (num(d.c) || 0) * unitC(d.cunit);
    const zc = () => capZ(w, C, num(d.esr) || 0, (num(d.esl) || 0) * 1e-9);
    const zl = () => indZ(w, L, num(d.dcr) || 0, (num(d.cpar) || 0) * 1e-12);
    switch (d.type) {
      case 'R': ck.z(inN, outN, cx(R, 0)); break;
      case 'C': ck.short(inN, outN); ck.z(outN, refNode, zc()); break;
      case 'RC': ck.z(inN, outN, cx(R, 0)); ck.z(outN, refNode, zc()); break;
      case 'LC': ck.z(inN, outN, zl()); ck.z(outN, refNode, zc()); break;
      case 'CL': ck.z(inN, refNode, zc()); ck.z(inN, outN, zl()); break;
      default: ck.short(inN, outN);
    }
  }

  /* ---------- Complétude de la configuration ---------- */
  function termComplete(w) {
    const e = w.cl || {};
    if (w.func === 'gnd') {
      if (e.ref !== 'chassis') return true;
      if (e.type === 'Z') return pos(e.z);
      if (e.type === 'C') return pos(e.c);
      return true;
    }
    if (w.func === 'shield') return e.type === 'zero' || pos(e.z);
    switch (e.type) {
      case 'open': return true;
      case 'R': return pos(e.r);
      case 'Z': return pos(e.z);
      case 'C': return pos(e.c);
      case 'RCs': case 'RCp': return pos(e.r) && pos(e.c);
      case 'LC': return pos(e.l) && pos(e.c);
      default: return false;
    }
  }
  function filterComplete(d) {
    if (!d) return false;
    switch (d.type) {
      case 'none': return true;
      case 'R': return pos(d.r);
      case 'C': return pos(d.c);
      case 'RC': return pos(d.r) && pos(d.c);
      case 'LC': case 'CL': return pos(d.l) && pos(d.c);
      default: return false;
    }
  }
  function loadComplete(d) {
    switch (d.loadType) {
      case 'R': case 'Z': return pos(d.loadR);
      case 'C': return pos(d.loadC);
      case 'L': return pos(d.loadL);
      case 'RCs': case 'RCp': return pos(d.loadR) && pos(d.loadC);
      case 'LC': return pos(d.loadL) && pos(d.loadC);
      default: return false;
    }
  }
  const isShieldedInternal = t => ['Fil blindé', 'Nappe blindée', 'Câble coaxial blindé'].includes(t);
  const isRibbon = t => t === 'Nappe non blindée' || t === 'Nappe blindée';

  // Faisceau blindé, modèle forfaitaire : une atténuation A (dB) ≥ 0 doit être saisie.
  function shieldComplete(H) { return H.shield !== 'shielded' || nonneg(H.shieldAttenDb); }
  function shieldFactor(H) { return H.shield === 'shielded' ? Math.pow(10, -num(H.shieldAttenDb) / 20) : 1; }
  function checkConfig(cfg) {
    const H = cfg.harness, I = cfg.internal, J = cfg.injection, N = H.wires.length;
    if (cfg.model === 'approx') {
      const blocksA = {
        injection: pos(J.iBCI) && nonneg(J.probePos) && num(J.probePos) <= num(H.length),
        harness: pos(H.length) && N >= 1 && (H.bundleMode === 'manual' ? pos(H.bundleD) : pos(H.isoD)) && (H.rMode === 'none' || (pos(H.wireD) && !(num(H.wireD) > num(H.isoD)))) && shieldComplete(H) && H.wires.every(termComplete),
        connector: cfg.connector.mode === 'transparent' || H.wires.every((w, k) => w.func === 'gnd' || w.func === 'shield' || filterComplete(cfg.connector.wires[k])),
        internal: I.present === false || (!!I.type && pos(I.length) && pos(I.height) && (isRibbon(I.type) ? pos(I.pitch) : (I.bundleMode !== 'manual' || pos(I.bundleD)))),
        cardConnector: true,
        filter: H.wires.every((w, k) => w.func === 'gnd' || w.func === 'shield' || (filterComplete(cfg.filters[k]) && loadComplete(cfg.filters[k]))) &&
          (cfg.cardGnd.type === 'direct' || cfg.cardGnd.type === 'open' || (cfg.cardGnd.type === 'Z' && pos(cfg.cardGnd.z)) || (cfg.cardGnd.type === 'C' && pos(cfg.cardGnd.c)))
      };
      const msgs = [];
      const hasReturn = H.wires.some(w => (w.func === 'gnd' && w.cl && w.cl.ref === 'chassis') || w.func === 'shield');
      if (!hasReturn) msgs.push('0 V de la charge non relié au châssis : à basse fréquence, le courant BCI ne revient que par la capacité du faisceau vers le plan. Imposer le courant y demande une f.é.m. irréaliste ; ces points sont signalés « inatteignable ».');
      if (H.shield === 'shielded') msgs.push(`Faisceau extérieur blindé : atténuation forfaitaire de ${num(H.shieldAttenDb)} dB appliquée aux fils intérieurs (approximation, voir les hypothèses).`);
      if (N > MAX_WIRES.approx) msgs.push(`Plus de ${MAX_WIRES.approx} conducteurs : le calcul approché est limité à ${MAX_WIRES.approx} conducteurs.`);
      return { blocks: blocksA, ready: Object.values(blocksA).every(Boolean) && N <= MAX_WIRES.approx, messages: msgs };
    }
    const blocks = {
      injection: pos(J.iBCI) && nonneg(J.probePos) && num(J.probePos) <= num(H.length) && (J.mode !== 'calibrated' || nonneg(J.zp)) && (J.limitDb === '' || J.limitDb === undefined || Number.isFinite(num(J.limitDb))),
      harness: pos(H.length) && (H.bundleMode !== 'manual' || pos(H.bundleD)) && pos(H.wireD) && pos(H.isoD) && num(H.isoD) >= num(H.wireD) && pos(H.epsr) && num(H.epsr) >= 1 && N >= 1 && shieldComplete(H) && H.wires.every(termComplete),
      connector: cfg.connector.mode === 'transparent' || H.wires.every((w, k) => w.func === 'gnd' || w.func === 'shield' || filterComplete(cfg.connector.wires[k])),
      internal: I.present === false || (!!I.type && pos(I.length) && pos(I.height) && (isRibbon(I.type) ? pos(I.pitch) : (I.bundleMode !== 'manual' || pos(I.bundleD)))),
      cardConnector: true,
      filter: H.wires.every((w, k) => w.func === 'gnd' || w.func === 'shield' || (filterComplete(cfg.filters[k]) && loadComplete(cfg.filters[k]))) &&
        (cfg.cardGnd.type === 'direct' || cfg.cardGnd.type === 'open' || (cfg.cardGnd.type === 'Z' && pos(cfg.cardGnd.z)) || (cfg.cardGnd.type === 'C' && pos(cfg.cardGnd.c)))
    };
    const messages = [];
    if (H.shield === 'shielded') messages.push(`Faisceau extérieur blindé : atténuation forfaitaire de ${num(H.shieldAttenDb)} dB appliquée aux fils intérieurs (approximation, voir les hypothèses).`);
    if (N > MAX_WIRES.detailed) messages.push(`Plus de ${MAX_WIRES.detailed} conducteurs : le modèle détaillé est limité à ${MAX_WIRES.detailed} conducteurs (temps de calcul). Utiliser le modèle approché, qui accepte jusqu'à ${MAX_WIRES.approx} conducteurs.`);
    else if (N > 25) messages.push(`${N} conducteurs en modèle détaillé : le calcul peut prendre quelques secondes.`);
    if (!H.wires.some(w => w.func === 'gnd')) messages.push('Aucun conducteur GND / 0 V : les tensions carte sont référencées au 0 V carte seul.');
    if (pos(H.length) && nonneg(J.probePos) && num(J.probePos) > num(H.length)) messages.push('La pince est positionnée au-delà de la longueur du harnais.');
    return { blocks, ready: Object.values(blocks).every(Boolean) && N <= MAX_WIRES.detailed, messages };
  }

  /* ---------- Diamètres effectifs (mm) ---------- */
  function harnessBundleD(H) { return H.bundleMode === 'manual' ? num(H.bundleD) : autoBundleD(H.wires.length, num(H.isoD)); }
  function internalBundleD(I, N, isoD) { return I.bundleMode === 'manual' ? num(I.bundleD) : autoBundleD(N, isoD); }

  /* ---------- Préparation des lignes ---------- */
  function prepare(cfg) {
    const H = cfg.harness, I = cfg.internal, N = H.wires.length;
    if (cfg.model === 'approx') {
      // Mode approché : faisceau équivalent (cylindre de diamètre D à la hauteur h du plan), ligne de mode commun
      // L′cm = μ0/2π·acosh(h/R), C′cm = 2πε0/acosh(h/R), répartie également entre les N conducteurs
      // (lignes non couplées : L′w = N·L′cm, C′w = C′cm/N, v = c). Couplages entre conducteurs négligés.
      const mk = (n, Dmm, h, label) => {
        const R = Dmm / 2000;
        if (!(R > 0)) return { errors: [`${label} : diamètre du faisceau non défini.`], warnings: [] };
        if (!(h > R * 1.0001)) return { errors: [`${label} : la hauteur au-dessus du plan doit être supérieure au rayon du faisceau.`], warnings: [] };
        const a = Math.acosh(h / R), Lcm = MU0 / (2 * Math.PI) * a, Ccm = 2 * Math.PI * EPS0 / a;
        const Lw = n * Lcm, Cw = Ccm / n, v = C0, y = 1 / (v * Lw);
        const Yc = Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => (i === j ? y : 0)));
        const eye = Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => (i === j ? 1 : 0)));
        return { Yc, v, n, Lcm, Ccm, Zcm: v * Lcm, Lw, Cw, Zw: v * Lw, D: Dmm, h, eig: { vals: Array(n).fill(Cw), T: eye }, errors: [], warnings: [] };
      };
      const PH = mk(N, harnessBundleD(H), 0.05, 'Harnais extérieur');
      let PI = null;
      if (I.present !== false && pos(I.length)) {
        const Dint = isRibbon(I.type) ? N * num(I.pitch) / 2 : internalBundleD(I, N, num(H.isoD)); // nappe : rayon équivalent = largeur / 4
        PI = mk(N, Dint, num(I.height) / 1000, 'Câblage interne');
      }
      return { PH, PI, approx: true, rw: num(H.wireD) / 2000, errors: [...PH.errors, ...((PI && PI.errors) || [])], warnings: [] };
    }
    const rw = num(H.wireD) / 2000, iso = num(H.isoD) / 1000;
    const geoH = bundlePositions(N, harnessBundleD(H) / 1000, iso, 0.05);
    const PH = lineParams(geoH, rw, num(H.epsr) || 1);
    let PI = null;
    if (I.present !== false && pos(I.length)) {
      const h = num(I.height) / 1000;
      const geoI = isRibbon(I.type) ? rowPositions(N, num(I.pitch) / 1000, h) : bundlePositions(N, internalBundleD(I, N, num(H.isoD)) / 1000, iso, h);
      PI = lineParams(geoI, rw, 1);
    }
    return { PH, PI, rw, errors: [...(PH.errors || []), ...((PI && PI.errors) || [])], warnings: [...(PH.warnings || []), ...((PI && PI.warnings) || [])] };
  }

  /* ---------- Construction et résolution à une fréquence ---------- */
  function buildCircuit(cfg, prep, f) {
    const w = 2 * Math.PI * f, H = cfg.harness, N = H.wires.length, ck = new Circuit();
    const Rp = H.rMode === 'none' ? 0 : rPrime(f, prep.rw);
    const L0 = ck.node(), lw = ck.nodes(N);
    const gndIdx = H.wires.findIndex(x => x.func === 'gnd');
    // Terminaisons côté charge
    H.wires.forEach((wv, k) => {
      const e = wv.cl || {};
      if (wv.func === 'gnd') {
        ck.short(lw[k], L0);
        if (e.ref === 'chassis') {
          if (e.type === 'direct') ck.short(L0, 0);
          else if (e.type === 'Z') ck.z(L0, 0, cx(num(e.z), 0));
          else if (e.type === 'C') ck.z(L0, 0, capZ(w, num(e.c) * unitC(e.cunit)));
        }
      } else if (wv.func === 'shield') {
        if (e.type === 'zero') ck.short(lw[k], 0); else ck.z(lw[k], 0, cx(num(e.z), 0));
      } else ck.z(lw[k], L0, termZ(e, w));
    });
    // Harnais : charge → pince → connecteur équipement
    const Lh = num(H.length), xp = Math.min(Math.max(num(cfg.injection.probePos), 0), Lh);
    let pa = lw;
    if (Lh - xp > 1e-9) { pa = ck.nodes(N); ck.line(lw, pa, Lh - xp, prep.PH, w, Rp); }
    const pb = ck.nodes(N);
    for (let k = 0; k < N; k++) ck.vsrc(pa[k], pb[k]);
    let ce = pb;
    if (xp > 1e-9) { ce = ck.nodes(N); ck.line(pb, ce, xp, prep.PH, w, Rp); }
    // Connecteur entrée équipement
    const ci = ck.nodes(N);
    const refEq = gndIdx >= 0 ? ci[gndIdx] : 0;
    H.wires.forEach((wv, k) => {
      if (wv.func === 'shield') { ck.short(ce[k], 0); ck.short(ci[k], 0); return; } // reprise blindage au châssis équipement
      const d = cfg.connector.wires[k];
      if (cfg.connector.mode === 'transparent' || wv.func === 'gnd' || !d) ck.short(ce[k], ci[k]);
      else stampFilter(ck, ce[k], ci[k], d, d.cref === 'chassis' ? 0 : refEq, w);
    });
    // Câblage interne
    let cd = ci;
    if (prep.PI) { cd = ck.nodes(N); ck.line(ci, cd, num(cfg.internal.length), prep.PI, w, Rp); }
    // Carte
    const C0n = ck.node(), out = Array(N).fill(null);
    H.wires.forEach((wv, k) => {
      if (wv.func === 'gnd') { ck.short(cd[k], C0n); return; }
      if (wv.func === 'shield') return;
      const d = cfg.filters[k], cl = ck.node();
      stampFilter(ck, cd[k], cl, d, d.cref === 'chassis' ? 0 : C0n, w);
      ck.z(cl, C0n, loadZ(d, w));
      out[k] = cl;
    });
    const g = cfg.cardGnd;
    if (g.type === 'direct') ck.short(C0n, 0);
    else if (g.type === 'Z') ck.z(C0n, 0, cx(num(g.z), 0));
    else if (g.type === 'C') ck.z(C0n, 0, capZ(w, num(g.c) * unitC(g.cunit)));
    return { ck, nodes: { L0, lw, pa, pb, ce, ci, cd, C0n, out }, gndIdx };
  }

  const RJIG = 100; // boucle de calibration 2 × 50 Ω (hypothèse déclarée)

  function point(cfg, prep, f) {
    const { ck, nodes, gndIdx } = buildCircuit(cfg, prep, f);
    const J = cfg.injection, Ibci = num(J.iBCI) / 1000, N = cfg.harness.wires.length;
    let sol, factor = 1, E, Zp;
    const calibrated = J.mode === 'calibrated' && !prep.approx;
    if (calibrated) {
      Zp = cx(num(J.zp) || 0, 0); E = cx(Ibci * (RJIG + Zp.r), 0);
      sol = ck.solve(E, Zp); if (!sol) return null;
    } else {
      Zp = cx(0, 0); E = cx(1, 0);
      sol = ck.solve(E, Zp); if (!sol) return null;
    }
    let Itot = cx(0, 0); for (let k = 0; k < N; k++) Itot = add(Itot, sol.I(k));
    const iMag = abs(Itot);
    let limited = false;
    if (calibrated) {
      if (J.limitDb !== '' && J.limitDb !== undefined && Number.isFinite(num(J.limitDb))) {
        const Ilim = Ibci * Math.pow(10, num(J.limitDb) / 20);
        if (iMag > Ilim) { factor = Ilim / iMag; limited = true; }
      }
    } else {
      if (!(iMag > 1e-15)) return null; // aucun chemin de retour pour le courant de mode commun
      factor = Ibci / iMag;
    }
    // Faisceau blindé (modèle forfaitaire) : le courant BCI circule sur le blindage ; tout ce qui atteint les fils
    // intérieurs est réduit du facteur ks = 10^(−A/20), A étant l'atténuation saisie.
    const shielded = cfg.harness.shield === 'shielded', ks = shieldFactor(cfg.harness);
    const V = n => scl(sol.V(n), factor * ks);
    const vC0 = V(nodes.C0n);
    const wires = cfg.harness.wires.map((wv, k) => {
      const r = { func: wv.func, iProbe: abs(sol.I(k)) * factor * ks, iShare: abs(sol.I(k)) / iMag * ks, vConnector: abs(V(nodes.ce[k])), vCardIn: abs(V(nodes.cd[k])) };
      r.vConnectorDiff = gndIdx >= 0 ? abs(sub(V(nodes.ce[k]), V(nodes.ce[gndIdx]))) : NaN;
      r.vCardInDiff = gndIdx >= 0 ? abs(sub(V(nodes.cd[k]), V(nodes.cd[gndIdx]))) : NaN;
      if (nodes.out[k] != null) r.vout = abs(sub(V(nodes.out[k]), vC0));
      return r;
    });
    const Ieff = iMag * factor;
    return { f, Ibundle: Ieff, E: shielded ? NaN : (calibrated ? abs(E) * factor : factor), limited, shielded, ks, vCardGnd: abs(vC0), wires, dim: sol.dim };
  }

  // Mode approché : fréquence au-delà de laquelle la longueur dépasse λ/10 (propagation non représentée).
  function approxLimit(cfg) { const L = num(cfg.harness.length); return L > 0 ? C0 / (10 * L) : Infinity; }
  function bandLimits(band) { return band === 'cs114' ? [10e3, 200e6] : [10e3, 400e6]; }
  function frequencies(band, n = 601) {
    const [f0, f1] = bandLimits(band);
    return Array.from({ length: n }, (_, k) => f0 * Math.pow(f1 / f0, k / (n - 1)));
  }
  function sweep(cfg, n = 601) {
    const prep = prepare(cfg);
    if (prep.errors.length) return { errors: prep.errors, warnings: prep.warnings, points: [] };
    const points = frequencies(cfg.injection.band, n).map(f => point(cfg, prep, f)).filter(Boolean);
    return { errors: [], warnings: prep.warnings, points, prep };
  }

  /* ---------- Auto-tests analytiques ---------- */
  function selfTests() {
    const out = [];
    const rel = (a, b) => Math.abs(a - b) / Math.max(Math.abs(b), 1e-300);
    // T1 : Z0 d'un fil au-dessus du plan vs formule exacte 60·acosh(h/r)
    {
      const h = 0.05, rw = 0.5e-3, P = lineParams({ pts: [{ x: 0, y: h }], warnings: [] }, rw, 1);
      const z0 = 1 / P.Yc[0][0], ref = 60 * Math.acosh(h / rw);
      out.push({ name: 'Impédance caractéristique fil/plan (h = 5 cm, r = 0,5 mm)', ok: rel(z0, ref) < 1e-3, detail: `Z0 = ${z0.toFixed(2)} Ω ; 60·acosh(h/r) = ${ref.toFixed(2)} Ω` });
    }
    // T2 : ligne court-circuitée excitée en série : I = E / (j Z0 tan βl)
    {
      const h = 0.05, rw = 0.5e-3, P = lineParams({ pts: [{ x: 0, y: h }], warnings: [] }, rw, 1), Z0 = 1 / P.Yc[0][0], l = 1;
      let worst = 0;
      for (const f of [1e5, 10e6, 60e6, 120e6, 333e6]) {
        const w = 2 * Math.PI * f, ck = new Circuit(), n1 = ck.node(), n2 = ck.node();
        ck.vsrc(0, n1); ck.line([n1], [n2], l, P, w); ck.short(n2, 0);
        const s = ck.solve(cx(1, 0), cx(0, 0)), I = abs(s.I(0)), ref = 1 / Math.abs(Z0 * Math.tan(w * l / C0));
        worst = Math.max(worst, rel(I, ref));
      }
      out.push({ name: 'Ligne court-circuitée : courant vs 1/|Z0·tan βl| (0,1 à 333 MHz)', ok: worst < 1e-6, detail: `écart relatif max = ${worst.toExponential(2)}` });
    }
    // T3 : milieu homogène : L′C′ = I/v²
    {
      const geo = ringPositions(3, 0.01, 1.4e-3, 0.05), P = lineParams(geo, 0.4e-3, 1);
      let err = 0;
      for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) { let s = 0; for (let k = 0; k < 3; k++) s += P.L[i][k] * P.C[k][j]; err = Math.max(err, Math.abs(s * C0 * C0 - (i === j ? 1 : 0))); }
      out.push({ name: 'Cohérence L′·C′ = I / c² (faisceau 3 fils)', ok: err < 1e-9, detail: `écart max = ${err.toExponential(2)}` });
    }
    // T4 : symétrie — deux conducteurs identiques et terminaisons identiques : aucune conversion MC → MD
    {
      const cfg = {
        injection: { iBCI: 100, mode: 'current', probePos: 0.05, band: 'do160', zp: '', limitDb: '' },
        harness: { length: 1, bundleD: 10, bundleMode: 'manual', wireD: 0.8, isoD: 1.4, epsr: 1, shield: 'none', wires: [0, 1].map(() => ({ func: 'signal', cl: { type: 'R', r: 50 } })) },
        connector: { mode: 'transparent', wires: [] },
        internal: { type: 'Fil non blindé', length: 0.3, height: 10, bundleD: 5, bundleMode: 'manual', pitch: 1.27 },
        cardGnd: { type: 'direct' },
        filters: [0, 1].map(() => ({ type: 'C', c: 1, cunit: 'nF', cref: 'gnd', loadType: 'R', loadR: 1000 }))
      };
      const prep = prepare(cfg); let worst = 0;
      for (const f of [1e5, 30e6, 150e6]) {
        const { ck, nodes } = buildCircuit(cfg, prep, f), s = ck.solve(cx(1, 0), cx(0, 0));
        const v1 = s.V(nodes.out[0]), v2 = s.V(nodes.out[1]);
        worst = Math.max(worst, abs(sub(v1, v2)) / abs(v1));
      }
      out.push({ name: 'Symétrie : 2 fils identiques → tension différentielle nulle', ok: worst < 1e-9, detail: `|V1 − V2| / |V1| max = ${worst.toExponential(2)}` });
    }
    // T5 : basse fréquence — fil court-circuité aux deux bouts : I ≈ E / (jω L′ l)
    {
      const h = 0.05, rw = 0.5e-3, P = lineParams({ pts: [{ x: 0, y: h }], warnings: [] }, rw, 1), f = 10e3, w = 2 * Math.PI * f;
      const ck = new Circuit(), n1 = ck.node(), n2 = ck.node(); ck.vsrc(0, n1); ck.line([n1], [n2], 1, P, w); ck.short(n2, 0);
      const I = abs(ck.solve(cx(1, 0), cx(0, 0)).I(0)), ref = 1 / (w * P.L[0][0]);
      out.push({ name: 'Limite basse fréquence : boucle inductive (10 kHz)', ok: rel(I, ref) < 1e-4, detail: `I = ${(I * 1e3).toFixed(4)} mA/V ; 1/(ωL′l) = ${(ref * 1e3).toFixed(4)} mA/V` });
    }
    // T6 : réciprocité — même réseau, source quasi-idéale de courant (1 A via 1 MΩ) en A puis en B
    {
      const geo = ringPositions(3, 0.01, 1.4e-3, 0.05), P = lineParams(geo, 0.4e-3, 1), w = 2 * Math.PI * 47e6;
      const trans = injectAtA => {
        const ck = new Circuit(), a = ck.nodes(3), b = ck.nodes(3);
        ck.line(a, b, 1, P, w);
        a.forEach((n, k) => ck.z(n, 0, cx(50 + 10 * k, 0)));
        b.forEach((n, k) => ck.z(n, 0, cx(20 + 5 * k, 0)));
        const A = a[0], B = b[2], inj = injectAtA ? A : B, meas = injectAtA ? B : A, s = ck.node();
        ck.vsrc(0, s); ck.z(s, inj, cx(1e6, 0));
        ck.z(meas, 0, cx(1e6, 0)); // même charge 1 MΩ sur le nœud de mesure : réseaux identiques
        return ck.solve(cx(1e6, 0), cx(0, 0)).V(meas);
      };
      const vAB = trans(true), vBA = trans(false), e = abs(sub(vAB, vBA)) / abs(vAB);
      out.push({ name: 'Réciprocité du réseau MTL (transimpédance A→B = B→A, 47 MHz)', ok: e < 1e-9, detail: `écart relatif = ${e.toExponential(2)}` });
    }
    // T7 : mode approché — N lignes non couplées en parallèle = ligne de mode commun du faisceau équivalent
    {
      const cfg = { model: 'approx', harness: { length: 1, bundleMode: 'manual', bundleD: 10, isoD: 1.4, wires: [{}, {}, {}] }, internal: {} };
      const P = prepare(cfg).PH; let worst = 0;
      const Psingle = { Yc: [[1 / P.Zcm]], v: C0 };
      for (const f of [1e5, 40e6, 160e6]) {
        const w = 2 * Math.PI * f;
        const a = new Circuit(), n1 = a.nodes(3), n2 = a.nodes(3), s1 = a.node();
        a.vsrc(0, s1); n1.forEach(n => a.short(n, s1)); a.line(n1, n2, 1, P, w); n2.forEach(n => a.z(n, 0, cx(30, 0)));
        const b = new Circuit(), m1 = b.node(), m2 = b.node();
        b.vsrc(0, m1); b.line([m1], [m2], 1, Psingle, w); b.z(m2, 0, cx(10, 0));
        const Ia = abs(a.solve(cx(1, 0), cx(0, 0)).I(0)), Ib = abs(b.solve(cx(1, 0), cx(0, 0)).I(0));
        worst = Math.max(worst, rel(Ia, Ib));
      }
      out.push({ name: 'Mode approché : 3 conducteurs équivalents = ligne de mode commun du faisceau', ok: worst < 1e-9, detail: `écart relatif max = ${worst.toExponential(2)} ; Z0 mode commun = ${P.Zcm.toFixed(1)} Ω` });
    }
    // T8 : ligne avec pertes court-circuitée : I = E·Y0·coth(γl), Y0 = √(y/z), γ = √(zy), z = R′ + jωL′, y = jωC′
    {
      const h = 0.05, rw = 0.4e-3, P = lineParams({ pts: [{ x: 0, y: h }], warnings: [] }, rw, 1), l = 1;
      let worst = 0;
      for (const f of [1e4, 1e6, 75e6, 150e6, 400e6]) {
        const w = 2 * Math.PI * f, Rp = rPrime(f, rw);
        const ck = new Circuit(), n1 = ck.node(), n2 = ck.node(); ck.vsrc(0, n1); ck.line([n1], [n2], l, P, w, Rp); ck.short(n2, 0);
        const I = abs(ck.solve(cx(1, 0), cx(0, 0)).I(0));
        const z = cx(Rp, w * P.L[0][0]), y = cx(0, w * P.C[0][0]), g = scl(csqrt(mul(z, y)), l), Y0 = csqrt(div(y, z));
        const e2 = cexp(scl(g, -2)), ref = abs(mul(Y0, div(add(cx(1, 0), e2), sub(cx(1, 0), e2))));
        worst = Math.max(worst, rel(I, ref));
      }
      out.push({ name: 'Ligne avec pertes (effet de peau) court-circuitée : I = E·Y0·coth(γl), de 10 kHz à 400 MHz', ok: worst < 1e-6, detail: `écart relatif max = ${worst.toExponential(2)}` });
    }
    // T9 : continuité — pertes infimes sur un faisceau de 3 fils = ligne sans pertes
    {
      const geo = ringPositions(3, 0.01, 1.4e-3, 0.05), P = lineParams(geo, 0.4e-3, 1); let worst = 0;
      for (const f of [1e5, 33e6, 123e6]) {
        const w = 2 * Math.PI * f, run = Rp => {
          const ck = new Circuit(), a = ck.nodes(3), b = ck.nodes(3), s1 = ck.node();
          ck.vsrc(0, s1); ck.z(s1, a[0], cx(50, 0)); a.slice(1).forEach(n => ck.z(n, 0, cx(10, 0)));
          ck.line(a, b, 1, P, w, Rp); b.forEach((n, k) => ck.z(n, 0, cx(20 + 30 * k, 0)));
          return ck.solve(cx(1, 0), cx(0, 0)).V(b[2]);
        };
        const v0 = run(0), v1 = run(1e-9); worst = Math.max(worst, abs(sub(v0, v1)) / abs(v0));
      }
      out.push({ name: 'Continuité : pertes infimes = ligne sans pertes (faisceau de 3 fils, modes couplés)', ok: worst < 1e-6, detail: `écart relatif max = ${worst.toExponential(2)}` });
    }
    // T10 : résistance réelle (R + Ls) // Cp : impédance analytique
    {
      const R = 10, ls = 1.2, cp = 0.3; let worst = 0;
      for (const f of [1e6, 1e8, 4e8]) {
        const w = 2 * Math.PI * f, z = realR(R, ls, cp, w);
        const zs = cx(R, w * ls * 1e-9), ref = div(zs, add(cx(1, 0), mul(zs, cx(0, w * cp * 1e-12))));
        worst = Math.max(worst, abs(sub(z, ref)) / abs(ref));
      }
      out.push({ name: 'Charge résistive réelle : (R + jωLs) en parallèle avec Cp', ok: worst < 1e-12, detail: `écart relatif max = ${worst.toExponential(2)}` });
    }
    // T11 : faisceau blindé, modèle forfaitaire — tensions carte réduites exactement de A dB
    {
      const mk = sh => ({ model: 'approx', injection: { iBCI: 100, band: 'do160', probePos: 0.05 },
        harness: { length: 1, bundleMode: 'auto', isoD: 1.4, wireD: 0.8, shield: sh, shieldAttenDb: 26, wires: [{ func: 'signal', cl: { type: 'R', r: 100 } }, { func: 'gnd', cl: { ref: 'chassis', type: 'direct' } }] },
        connector: { mode: 'transparent', wires: [] }, internal: { present: false }, cardGnd: { type: 'direct' },
        filters: [{ type: 'none', loadType: 'R', loadR: 50 }, { type: 'none' }] });
      const a = mk('none'), b = mk('shielded'), pa = prepare(a), pb = prepare(b); let worst = 0;
      for (const f of [1e5, 3e7, 2e8]) { const va = point(a, pa, f).wires[0].vout, vb = point(b, pb, f).wires[0].vout; worst = Math.max(worst, Math.abs(20 * Math.log10(va / vb) - 26)); }
      out.push({ name: 'Faisceau blindé (forfaitaire) : tension carte réduite de l\'atténuation saisie (26 dB)', ok: worst < 1e-9, detail: `écart max = ${worst.toExponential(2)} dB` });
    }
    return out;
  }

  const api = { shieldFactor, realR, rPrime, rDC, skinDepth, jacobiEig, MAX_WIRES, C0, MU0, EPS0, cx, abs, ringPositions, bundlePositions, packingRatio, autoBundleD, harnessBundleD, internalBundleD, rowPositions, lineParams, inductanceMatrix, Circuit, prepare, buildCircuit, point, sweep, frequencies, bandLimits, approxLimit, checkConfig, selfTests, termComplete, filterComplete, loadComplete, isShieldedInternal, isRibbon, RJIG };
  if (typeof module !== 'undefined' && module.exports) module.exports = api; else global.BCIEngine = api;
})(typeof window !== 'undefined' ? window : globalThis);
