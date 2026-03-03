#[derive(Debug)]
struct Record<'a> {
    kind: &'a str,
    iteration: usize,
    evaluations: usize,
    f: f64,
    projected_gradient: f64,
    x: Vec<f64>,
    g: Vec<f64>,
    task: &'a str,
}

fn sha256(input: &[u8]) -> String {
    const INITIAL: [u32; 8] = [
        0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab,
        0x5be0cd19,
    ];
    const K: [u32; 64] = [
        0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4,
        0xab1c5ed5, 0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe,
        0x9bdc06a7, 0xc19bf174, 0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f,
        0x4a7484aa, 0x5cb0a9dc, 0x76f988da, 0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
        0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967, 0x27b70a85, 0x2e1b2138, 0x4d2c6dfc,
        0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85, 0xa2bfe8a1, 0xa81a664b,
        0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070, 0x19a4c116,
        0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
        0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7,
        0xc67178f2,
    ];

    let bit_len = (input.len() as u64) * 8;
    let mut padded = input.to_vec();
    padded.push(0x80);
    while padded.len() % 64 != 56 {
        padded.push(0);
    }
    padded.extend_from_slice(&bit_len.to_be_bytes());

    let mut hash = INITIAL;
    for block in padded.chunks_exact(64) {
        let mut words = [0u32; 64];
        for (word, bytes) in words.iter_mut().zip(block.chunks_exact(4)) {
            *word = u32::from_be_bytes(bytes.try_into().unwrap());
        }
        for i in 16..64 {
            let s0 = words[i - 15].rotate_right(7)
                ^ words[i - 15].rotate_right(18)
                ^ (words[i - 15] >> 3);
            let s1 = words[i - 2].rotate_right(17)
                ^ words[i - 2].rotate_right(19)
                ^ (words[i - 2] >> 10);
            words[i] = words[i - 16]
                .wrapping_add(s0)
                .wrapping_add(words[i - 7])
                .wrapping_add(s1);
        }

        let [mut a, mut b, mut c, mut d, mut e, mut f, mut g, mut h] = hash;
        for i in 0..64 {
            let sum1 = h
                .wrapping_add(e.rotate_right(6) ^ e.rotate_right(11) ^ e.rotate_right(25))
                .wrapping_add((e & f) ^ ((!e) & g))
                .wrapping_add(K[i])
                .wrapping_add(words[i]);
            let sum0 = (a.rotate_right(2) ^ a.rotate_right(13) ^ a.rotate_right(22))
                .wrapping_add((a & b) ^ (a & c) ^ (b & c));
            h = g;
            g = f;
            f = e;
            e = d.wrapping_add(sum1);
            d = c;
            c = b;
            b = a;
            a = sum1.wrapping_add(sum0);
        }
        for (state, value) in hash.iter_mut().zip([a, b, c, d, e, f, g, h]) {
            *state = state.wrapping_add(value);
        }
    }

    hash.iter().map(|word| format!("{word:08x}")).collect()
}

#[test]
fn fixture_bytes_are_canonical_and_pinned() {
    let fixtures: [(&str, &[u8], &str); 6] = [
        (
            "audit-mixed-rollover.csv",
            include_bytes!("reference/fixtures/audit-mixed-rollover.csv"),
            "a0e0c36fa02cedcc8140277cd056aa631c25603881ea0c5ddf12fa4b855bb95e",
        ),
        (
            "driver1.csv",
            include_bytes!("reference/fixtures/driver1.csv"),
            "b15ffe542100c7ed0440ea1a3f6d535e10d029bb31400a5c47648627982e983c",
        ),
        (
            "driver2.csv",
            include_bytes!("reference/fixtures/driver2.csv"),
            "a09051311c9c94c4e7b4d8d4a8b46e793a5a679de6c8e8f2c33ece65a4055c53",
        ),
        (
            "dcsrch-cases.csv",
            include_bytes!("reference/fixtures/dcsrch-cases.csv"),
            "2de020f68ad1607d9dd72ba6ceaca94e88e3121c27ec40ef7eb52559ad2462de",
        ),
        (
            "driver3-large-n.csv",
            include_bytes!("reference/fixtures/driver3-large-n.csv"),
            "789b9ae5cdddf49d4bc5b73f842b303043811b8418a79ddbc7cb4a9f3c45c8ce",
        ),
        (
            "edge-mixed-bounds.csv",
            include_bytes!("reference/fixtures/edge-mixed-bounds.csv"),
            "28d98e7e8ae4e796b96f374966a65e6910cc3778b47161c672978bbfa4f62895",
        ),
    ];

    for (name, bytes, expected) in fixtures {
        assert!(
            !bytes.contains(&b'\r'),
            "{name} must use canonical LF endings"
        );
        assert_eq!(sha256(bytes), expected, "{name} SHA-256 mismatch");
    }
}

fn parse_vector(field: &str) -> Vec<f64> {
    let field = field.trim_matches('"');
    if field.is_empty() {
        return Vec::new();
    }
    field
        .split('|')
        .map(|value| value.trim().parse().expect("valid full-precision float"))
        .collect()
}

fn parse_fixture(text: &str) -> Vec<Record<'_>> {
    let mut lines = text.lines();
    assert_eq!(
        lines.next(),
        Some("record,iteration,evaluations,f,projected_gradient,x,g,task")
    );
    lines
        .map(|line| {
            let fields = line.splitn(8, ',').collect::<Vec<_>>();
            assert_eq!(fields.len(), 8, "malformed fixture row: {line}");
            Record {
                kind: fields[0],
                iteration: fields[1].parse().unwrap(),
                evaluations: fields[2].parse().unwrap(),
                f: fields[3].trim().parse().unwrap(),
                projected_gradient: fields[4].trim().parse().unwrap(),
                x: parse_vector(fields[5]),
                g: parse_vector(fields[6]),
                task: fields[7],
            }
        })
        .collect()
}

fn assert_trajectory<'a>(
    text: &'a str,
    dimension: usize,
    iterations: usize,
    evaluations: usize,
    terminal_task: &str,
    accepted_vectors: bool,
) -> Vec<Record<'a>> {
    let records = parse_fixture(text);
    let accepted = records
        .iter()
        .filter(|record| record.kind == "accepted")
        .collect::<Vec<_>>();
    assert_eq!(accepted.len(), iterations);

    let mut previous_f = f64::INFINITY;
    for (index, record) in accepted.iter().enumerate() {
        assert_eq!(record.iteration, index + 1);
        assert!(record.evaluations >= record.iteration);
        assert!(record.f.is_finite());
        assert!(record.projected_gradient.is_finite());
        assert!(record.f <= previous_f);
        previous_f = record.f;
        if accepted_vectors {
            assert_eq!(record.x.len(), dimension);
            assert_eq!(record.g.len(), dimension);
        } else {
            assert!(record.x.is_empty());
            assert!(record.g.is_empty());
        }
    }

    let initial = records.first().unwrap();
    assert_eq!(initial.kind, "initial");
    assert_eq!(initial.evaluations, 1);
    assert_eq!(initial.x.len(), dimension);
    assert_eq!(initial.g.len(), dimension);

    let terminal = records.last().unwrap();
    assert_eq!(terminal.kind, "terminal");
    assert_eq!(terminal.iteration, iterations);
    assert_eq!(terminal.evaluations, evaluations);
    assert_eq!(terminal.x.len(), dimension);
    assert_eq!(terminal.g.len(), dimension);
    assert_eq!(terminal.task, terminal_task);
    records
}

#[test]
fn official_driver1_fixture_is_complete() {
    let records = assert_trajectory(
        include_str!("reference/fixtures/driver1.csv"),
        25,
        23,
        28,
        "CONVERGENCE: REL_REDUCTION_OF_F_<=_FACTR*EPSMCH",
        true,
    );
    assert!((records.last().unwrap().f - 1.083_490_083_430_061_5e-9).abs() < 1e-24);
}

#[test]
fn official_driver2_fixture_is_complete() {
    let records = assert_trajectory(
        include_str!("reference/fixtures/driver2.csv"),
        25,
        46,
        53,
        "STOP: THE PROJECTED GRADIENT IS SUFFICIENTLY SMALL",
        true,
    );
    assert!(records.last().unwrap().projected_gradient <= 1e-10);
}

#[test]
fn deterministic_large_n_driver3_fixture_is_complete() {
    let records = assert_trajectory(
        include_str!("reference/fixtures/driver3-large-n.csv"),
        1000,
        49,
        58,
        "STOP: THE PROJECTED GRADIENT IS SUFFICIENTLY SMALL",
        false,
    );
    assert!(records.last().unwrap().projected_gradient <= 1e-10);
}

#[test]
fn independent_mixed_bound_fixture_records_projection() {
    let records = assert_trajectory(
        include_str!("reference/fixtures/edge-mixed-bounds.csv"),
        4,
        2,
        3,
        "CONVERGENCE: NORM_OF_PROJECTED_GRADIENT_<=_PGTOL",
        true,
    );
    assert_eq!(records.first().unwrap().x, [0.0, 4.0, 7.0, 0.0]);
    let final_x = &records.last().unwrap().x;
    for (actual, expected) in final_x.iter().zip([1.0, 2.0, -2.0, 3.0]) {
        assert!((actual - expected).abs() < 1e-12);
    }
}

#[test]
fn audit_fixture_covers_all_bounds_ties_and_rollover() {
    let records = assert_trajectory(
        include_str!("reference/fixtures/audit-mixed-rollover.csv"),
        8,
        8,
        11,
        "CONVERGENCE: NORM_OF_PROJECTED_GRADIENT_<=_PGTOL",
        true,
    );
    assert_eq!(
        records.first().unwrap().x,
        [0.0, 8.0, 2.0, -6.0, 1.5, 0.0, 0.0, 0.0]
    );
    let final_x = &records.last().unwrap().x;
    assert_eq!(final_x[4], 1.5);
    assert_eq!(&final_x[5..], &[1.0, -0.5, 0.5]);
}

#[test]
fn official_dcsrch_fixture_records_terminal_branches() {
    let rows = include_str!("reference/fixtures/dcsrch-cases.csv")
        .lines()
        .collect::<Vec<_>>();
    assert_eq!(rows[0], "case,step,task");
    assert_eq!(rows.len(), 4);
    assert!(rows[1].ends_with("WARNING: STP = STPMAX"));
    assert!(rows[2].ends_with("WARNING: STP = STPMIN"));
    assert!(rows[3].ends_with("CONVERGENCE"));
}
