//! FND-3 §3.2 — STL import + the generalized winding number. Raw
//! ray-parity fails silently and catastrophically on non-watertight or
//! degenerate meshes (whole regions flipped); the winding number (Jacobson
//! 2013) degrades gracefully — inside/outside by thresholding at 0.5, no
//! mesh repair. This S9 build evaluates the winding number EXACTLY, per
//! triangle, via the van Oosterom & Strackee (1983) atan2 solid-angle
//! formula in fixed mesh order (the Barnes-Hut approximation of §3.2 is a
//! performance tier that rides a later wave; exactness is what the as-built
//! record pins). Triangle-touch gating uses the Akenine-Möller SAT test.

use super::Geom3dError;

/// A triangle soup in authoring (Cartesian) coordinates. `tris` keeps the
/// file's triangle order — the winding sum iterates it fixed (§3.5).
/// Degenerate (zero-area) triangles are DROPPED at parse with the count
/// recorded here — never a silent skip.
#[derive(Debug, Clone, PartialEq)]
pub struct TriMesh {
    pub tris: Vec<[[f64; 3]; 3]>,
    pub dropped_degenerate: usize,
}

impl TriMesh {
    /// Parse BOTH binary STL (80-byte header + u32 LE count + 50-byte
    /// records) and ASCII STL. Disambiguation prefers the binary
    /// record-count size identity `len == 84 + 50·count` — a binary file
    /// may legally begin with "solid" in its free-text header, so the
    /// leading token alone is NOT trusted (the classic importer trap).
    pub fn from_stl_bytes(bytes: &[u8]) -> Result<TriMesh, Geom3dError> {
        if bytes.len() >= 84 {
            let count = u32::from_le_bytes([bytes[80], bytes[81], bytes[82], bytes[83]]) as u64;
            if 84 + 50 * count == bytes.len() as u64 {
                return Self::parse_binary(bytes, count as usize);
            }
        }
        let text = std::str::from_utf8(bytes)
            .map_err(|e| Geom3dError::StlParse(format!("not binary-sized and not UTF-8: {e}")))?;
        if text.trim_start().starts_with("solid") {
            Self::parse_ascii(text)
        } else {
            Err(Geom3dError::StlParse(
                "neither the binary size identity (84 + 50·count) nor an ASCII 'solid' header"
                    .to_string(),
            ))
        }
    }

    fn parse_binary(bytes: &[u8], count: usize) -> Result<TriMesh, Geom3dError> {
        let mut mesh = TriMesh {
            tris: Vec::with_capacity(count),
            dropped_degenerate: 0,
        };
        for k in 0..count {
            // Record: 12 bytes normal (ignored — recomputable), 3×12 bytes
            // vertices, 2 bytes attribute.
            let base = 84 + 50 * k + 12;
            let mut tri = [[0.0f64; 3]; 3];
            for (v, vert) in tri.iter_mut().enumerate() {
                for (c, coord) in vert.iter_mut().enumerate() {
                    let off = base + 12 * v + 4 * c;
                    let raw = f32::from_le_bytes([
                        bytes[off],
                        bytes[off + 1],
                        bytes[off + 2],
                        bytes[off + 3],
                    ]);
                    *coord = f64::from(raw);
                }
            }
            mesh.push_tri(tri)?;
        }
        Ok(mesh)
    }

    fn parse_ascii(text: &str) -> Result<TriMesh, Geom3dError> {
        let mut mesh = TriMesh {
            tris: Vec::new(),
            dropped_degenerate: 0,
        };
        let mut verts: Vec<[f64; 3]> = Vec::new();
        let mut toks = text.split_whitespace();
        while let Some(tok) = toks.next() {
            if tok != "vertex" {
                continue;
            }
            let mut v = [0.0f64; 3];
            for coord in &mut v {
                let t = toks.next().ok_or_else(|| {
                    Geom3dError::StlParse("truncated vertex in ASCII STL".to_string())
                })?;
                *coord = t.parse::<f64>().map_err(|e| {
                    Geom3dError::StlParse(format!("bad vertex coordinate {t:?}: {e}"))
                })?;
            }
            verts.push(v);
        }
        if verts.is_empty() || !verts.len().is_multiple_of(3) {
            return Err(Geom3dError::StlParse(format!(
                "ASCII STL vertex count {} is not a positive multiple of 3",
                verts.len()
            )));
        }
        for tri in verts.chunks_exact(3) {
            mesh.push_tri([tri[0], tri[1], tri[2]])?;
        }
        Ok(mesh)
    }

    /// Append with the degeneracy gate: non-finite vertices refuse loudly;
    /// exactly-zero-area triangles are dropped and counted.
    fn push_tri(&mut self, tri: [[f64; 3]; 3]) -> Result<(), Geom3dError> {
        for v in &tri {
            if !(v[0].is_finite() && v[1].is_finite() && v[2].is_finite()) {
                return Err(Geom3dError::StlParse(format!("non-finite vertex {v:?}")));
            }
        }
        let n = cross(sub(tri[1], tri[0]), sub(tri[2], tri[0]));
        if dot(n, n) == 0.0 {
            self.dropped_degenerate += 1;
        } else {
            self.tris.push(tri);
        }
        Ok(())
    }

    /// Non-empty + finite (meshes built by hand bypass the parser gates).
    pub fn validate(&self) -> Result<(), Geom3dError> {
        if self.tris.is_empty() {
            return Err(Geom3dError::BadSolid("empty triangle mesh".to_string()));
        }
        for (i, t) in self.tris.iter().enumerate() {
            for v in t {
                if !(v[0].is_finite() && v[1].is_finite() && v[2].is_finite()) {
                    return Err(Geom3dError::BadSolid(format!(
                        "mesh triangle {i} has non-finite vertex {v:?}"
                    )));
                }
            }
        }
        Ok(())
    }

    /// Area-weighted mean outward normal of the triangles touching an
    /// axis-aligned box — the mesh path's PLIC normal operand. The winding
    /// number of a watertight mesh is piecewise CONSTANT (0/1), so a
    /// central difference of the classifier field is analytically zero
    /// away from the surface and cannot orient a PLIC plane; the mesh's
    /// own facets carry the orientation instead: Σ (b−a)×(c−a) over
    /// touching triangles (each cross product = 2·area·n̂_outward — the
    /// weighting is free), fixed mesh order (§3.5). `None` when no
    /// triangle touches or the sum cancels to zero — the caller refuses.
    pub fn interface_normal_in_box(&self, bmin: [f64; 3], bmax: [f64; 3]) -> Option<[f64; 3]> {
        let mut acc = [0.0f64; 3];
        let mut touched = false;
        for t in &self.tris {
            if tri_touches_box(t, bmin, bmax) {
                touched = true;
                let n = cross(sub(t[1], t[0]), sub(t[2], t[0]));
                acc[0] += n[0];
                acc[1] += n[1];
                acc[2] += n[2];
            }
        }
        if touched && dot(acc, acc) > 0.0 {
            Some(acc)
        } else {
            None
        }
    }

    /// The EXACT generalized winding number at `p`: the per-triangle
    /// signed solid angle (van Oosterom & Strackee 1983),
    /// `Ω = 2·atan2(a·(b×c), |a||b||c| + (a·b)|c| + (b·c)|a| + (c·a)|b|)`
    /// with a, b, c the vertices relative to `p`, summed in fixed mesh
    /// order (§3.5) and divided by 4π. Integer-valued for watertight
    /// outward-oriented meshes (1 inside, 0 outside); degrades gracefully
    /// for imperfect ones. Inside ⇔ w > 0.5 (the FND-3 §3.2 threshold).
    pub fn winding_number(&self, p: [f64; 3]) -> f64 {
        let mut sum = 0.0f64;
        for t in &self.tris {
            let a = sub(t[0], p);
            let b = sub(t[1], p);
            let c = sub(t[2], p);
            let la = dot(a, a).sqrt();
            let lb = dot(b, b).sqrt();
            let lc = dot(c, c).sqrt();
            let num = dot(a, cross(b, c));
            let den = la * lb * lc + dot(a, b) * lc + dot(b, c) * la + dot(c, a) * lb;
            sum += 2.0 * num.atan2(den);
        }
        sum / (4.0 * std::f64::consts::PI)
    }
}

/// Akenine-Möller triangle vs axis-aligned-box SAT test (the §3.2
/// triangle-touch gate). Separation uses strict inequalities, so a
/// touching triangle counts as intersecting — conservative for the
/// pure-cell short-circuit (touching ⇒ not provably pure ⇒ sample).
pub fn tri_touches_box(tri: &[[f64; 3]; 3], bmin: [f64; 3], bmax: [f64; 3]) -> bool {
    let c = [
        0.5 * (bmin[0] + bmax[0]),
        0.5 * (bmin[1] + bmax[1]),
        0.5 * (bmin[2] + bmax[2]),
    ];
    let e = [
        0.5 * (bmax[0] - bmin[0]),
        0.5 * (bmax[1] - bmin[1]),
        0.5 * (bmax[2] - bmin[2]),
    ];
    let v = [sub(tri[0], c), sub(tri[1], c), sub(tri[2], c)];
    let f = [sub(v[1], v[0]), sub(v[2], v[1]), sub(v[0], v[2])];

    // 9 cross axes e_i × f_j (a zero cross axis projects everything to 0
    // and reports no separation — correct: it is a degenerate axis).
    for i in 0..3 {
        for fj in &f {
            let mut axis = [0.0f64; 3];
            let mut unit = [0.0f64; 3];
            unit[i] = 1.0;
            axis[0] = unit[1] * fj[2] - unit[2] * fj[1];
            axis[1] = unit[2] * fj[0] - unit[0] * fj[2];
            axis[2] = unit[0] * fj[1] - unit[1] * fj[0];
            let p0 = dot(axis, v[0]);
            let p1 = dot(axis, v[1]);
            let p2 = dot(axis, v[2]);
            let r = e[0] * axis[0].abs() + e[1] * axis[1].abs() + e[2] * axis[2].abs();
            let pmin = p0.min(p1).min(p2);
            let pmax = p0.max(p1).max(p2);
            if pmin > r || pmax < -r {
                return false;
            }
        }
    }
    // 3 box face axes.
    for i in 0..3 {
        let pmin = v[0][i].min(v[1][i]).min(v[2][i]);
        let pmax = v[0][i].max(v[1][i]).max(v[2][i]);
        if pmin > e[i] || pmax < -e[i] {
            return false;
        }
    }
    // Triangle plane vs box.
    let n = cross(f[0], f[1]);
    let d = dot(n, v[0]);
    let r = e[0] * n[0].abs() + e[1] * n[1].abs() + e[2] * n[2].abs();
    if d > r || d < -r {
        return false;
    }
    true
}

#[inline]
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
