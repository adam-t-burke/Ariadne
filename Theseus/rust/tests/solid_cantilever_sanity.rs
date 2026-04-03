use ndarray::Array2;
use theseus::solid_solve::solve_solid;
use theseus::solid_types::{
    SolidCache, SolidDofMap, SolidElementProps, SolidLoad, SolidMaterial, SolidProblem, SolidSupport,
};

#[test]
fn skewed_tet_balances_load_and_deflects_downward() {
    let node_positions = Array2::from_shape_vec(
        (4, 3),
        vec![
            0.0, 0.0, 0.0,
            2.0, 0.0, 0.0,
            1.0, 3.0, 0.0,
            0.0, 0.0, 4.0,
        ],
    )
    .unwrap();

    let supports = vec![
        SolidSupport { node_idx: 0, fixed_dofs: [true, true, true] },
        SolidSupport { node_idx: 1, fixed_dofs: [true, true, true] },
        SolidSupport { node_idx: 2, fixed_dofs: [true, true, true] },
    ];

    let loads = vec![SolidLoad {
        node_idx: 3,
        force: [0.0, 0.0, -1_000.0],
    }];

    let problem = SolidProblem {
        num_nodes: 4,
        num_elements: 1,
        materials: vec![SolidMaterial {
            e: 210e9,
            nu: 0.3,
            density: 7850.0,
            yield_stress: 250e6,
        }],
        element_props: vec![SolidElementProps { material_idx: 0 }],
        supports: supports.clone(),
        loads,
        node_positions,
        elements: vec![vec![0, 1, 2, 3]],
        dof_map: SolidDofMap::from_supports(4, &supports),
        include_self_weight: false,
        gravity: [0.0, 0.0, -9.81],
    };

    let mut cache = SolidCache::new(&problem).unwrap();
    let result = solve_solid(&mut cache, &problem).unwrap();

    let rz: f64 = (0..4).map(|n| result.reactions[n * 3 + 2]).sum();
    assert!((rz - 1_000.0).abs() < 1e-6, "reactions should balance applied load");

    let ux = result.displacements[9];
    let uy = result.displacements[10];
    let uz = result.displacements[11];
    assert!(ux.is_finite() && uy.is_finite() && uz.is_finite());
    assert!(uz < 0.0, "loaded node should move downward");

    for gp_vm in result.von_mises {
        for vm in gp_vm {
            assert!(vm.is_finite());
        }
    }
}
