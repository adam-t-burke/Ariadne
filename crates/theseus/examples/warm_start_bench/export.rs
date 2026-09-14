//! Export of the synthetic suite nets as JSON cases (`export_cases <out_dir>`).
//!
//! Writes one file per suite net in the same schema as the external cases
//! (`bench/external/cases/*.json`, see `external.rs`): `nodes` / `target`
//! at the funicular, `loads` per node, `q_ref`, `signs`, the loose box in
//! `bounds`, plus two extra fields the Grasshopper loader understands:
//! `bounds_snug` (the ×/÷1 box around `q_true`) and `target_variants`
//! (`jit2pctd`, `bump10pctd`: the perturbed targets the suite and the
//! presentation figures use). Because `target == nodes`, the files can also
//! be replayed through `warm_start_bench external`.

use crate::nets::{suite_nets, Net};
use ndarray::Array2;
use serde_json::{json, Value};
use std::path::Path;

/// Perturbed targets exported per net: (key, jitter, bump), see [`Net::target`].
const VARIANTS: &[(&str, f64, f64)] = &[("jit2pctd", 0.02, 0.0), ("bump10pctd", 0.0, 0.10)];

fn array_json(x: &Array2<f64>) -> Value {
    Value::Array(
        x.rows()
            .into_iter()
            .map(|r| json!([r[0], r[1], r[2]]))
            .collect(),
    )
}

fn bounds_json(lo: &[f64], hi: &[f64]) -> Value {
    json!({ "lo": lo, "hi": hi })
}

/// Nodal load on every node (zeros at the supports), row-major `n × 3`.
fn full_loads(net: &Net) -> Vec<[f64; 3]> {
    let mut loads = vec![[0.0; 3]; net.n_nodes()];
    for (i, &node) in net.free.iter().enumerate() {
        loads[node] = match &net.load_xyz {
            Some(v) => v[i],
            None => [0.0, 0.0, -net.loads[i]],
        };
    }
    loads
}

fn case_json(net: &Net) -> Option<Value> {
    let funicular = net.try_funicular()?;
    let nodes = net.full_positions(&funicular);
    let (loose_lo, loose_hi) = net.box_from_true(100.0, 100.0);
    let (snug_lo, snug_hi) = net.box_from_true(1.0, 1.0);
    let signs: Vec<i32> = net.q_true.iter().map(|q| q.signum() as i32).collect();
    let tension = signs.iter().filter(|&&s| s > 0).count();
    let compression = signs.len() - tension;

    let mut variants = serde_json::Map::new();
    for &(key, jitter, bump) in VARIANTS {
        let target = net.full_positions(&net.target(jitter, bump));
        variants.insert(key.to_string(), array_json(&target));
    }

    Some(json!({
        "name": net.name,
        "source": "warm_start_bench",
        "source_ref": "crates/theseus/examples/warm_start_bench/nets.rs",
        "description": format!(
            "Synthetic suite net '{}' at its funicular ({} tension / {} compression members).",
            net.name, tension, compression
        ),
        "nodes": array_json(&nodes),
        "edges": net.edges.iter().map(|&(a, b)| [a, b]).collect::<Vec<_>>(),
        "fixed": net.fixed,
        "loads": full_loads(net),
        "q_ref": net.q_true,
        "signs": signs,
        "target": array_json(&nodes),
        "bounds": bounds_json(&loose_lo, &loose_hi),
        "bounds_snug": bounds_json(&snug_lo, &snug_hi),
        "target_variants": Value::Object(variants),
        "extent": net.extent,
        "tie_edges": net.tie_edges,
        "sign_mix": { "tension": tension, "compression": compression },
    }))
}

pub fn cmd_export_cases(out_dir: &str) {
    let out = Path::new(out_dir);
    std::fs::create_dir_all(out).expect("create output directory");
    for net in suite_nets() {
        let Some(value) = case_json(&net) else {
            eprintln!("{}: funicular forward solve failed, skipped", net.name);
            continue;
        };
        let path = out.join(format!("{}.json", net.name));
        let text = serde_json::to_string(&value).expect("serialise");
        std::fs::write(&path, text).unwrap_or_else(|e| panic!("write {}: {e}", path.display()));
        println!(
            "wrote {} ({} nodes, {} edges, {} fixed)",
            path.display(),
            net.n_nodes(),
            net.edges.len(),
            net.fixed.len()
        );
    }
}
