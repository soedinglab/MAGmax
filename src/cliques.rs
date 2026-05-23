use std::collections::{HashMap, HashSet};
use petgraph::graph::{Graph};
use petgraph::graph::NodeIndex;
use petgraph::{Undirected};
use rayon::prelude::*;
use crate::assess::BinQuality;
use crate::merge::AniData;

// Get clique clusters
pub fn split_component_into_cliques(
    component: HashSet<u32>,
    ani_map: &HashMap<(u32, u32), AniData>,
    ani_cutoff: f32,
    aligned_frac: f32,
    id_to_name: &[String],
    bin_qualities: &HashMap<String, BinQuality>,
    isolate_genomes: &HashSet<String>,
    no_reassembly: bool,
) -> Vec<HashSet<u32>> {

    let adj = build_adj(
        &component,
        ani_map,
        ani_cutoff,
        aligned_frac,
    );

    let mut remaining = component;
    let mut subclusters: Vec<HashSet<u32>> = Vec::new();
    let mut cores: Vec<HashSet<u32>> = Vec::new();

    const MAX_CLIQUE_CORE_SIZE: usize = 80;

    // Split large components into cores to improve performance
    // It might miss some cliques spanning multiple cores
    // but since we compare final clusters later, it should 
    // not miss genomic pairs > ANI threshold and not influence the final results
    while remaining.len() > MAX_CLIQUE_CORE_SIZE {

        let core = 
            extract_core_chunk(&mut remaining, &adj, MAX_CLIQUE_CORE_SIZE);
        if core.len() <= 1 {
            // treat as singleton(s)
            subclusters.extend(core.into_iter().map(|id| HashSet::from([id])));
            continue;

        }
        cores.push(core);
    }
    if !cores.is_empty() {
        let all_cliques: Vec<HashSet<u32>> = cores
            .par_iter()
            .flat_map(|core| {
                let subgraph = build_subgraph_for_ids(core, &adj);

                let mut r: Vec<NodeIndex> = Vec::new();
                let mut p: Vec<NodeIndex> = subgraph.node_indices().collect();
                let mut x: Vec<NodeIndex> = Vec::new();
                let mut cliques_vec: Vec<Vec<u32>> = Vec::new();

                bron_kerbosch(&subgraph, &mut r, &mut p, &mut x, &mut cliques_vec);

                cliques_vec
                    .into_iter()
                    .map(|c| c.into_iter().collect::<HashSet<u32>>())
                    .collect::<Vec<_>>() // returned from this closure
            })
            .collect();
        
        let mut sorted_cliques: Vec<Vec<u32>> = all_cliques
            .into_iter()
            .map(|set| {
                let mut v: Vec<u32> = set.into_iter().collect();
                v.sort_unstable();
                v
            })
            .collect();

        sorted_cliques.sort_unstable();
        subclusters.extend(
            sorted_cliques
            .into_iter()
            .map(|c| c.into_iter().collect::<HashSet<u32>>()),
        );
    }

    if !remaining.is_empty() {
        if remaining.len() == 1 {
            subclusters.push(HashSet::from([*remaining.iter().next().unwrap()]));
        } else {
            let subgraph = 
                build_subgraph_for_ids(&remaining, &adj);

            let mut r: Vec<NodeIndex> = Vec::new();
            let mut p: Vec<NodeIndex> = subgraph.node_indices().collect();
            let mut x: Vec<NodeIndex> = Vec::new();
            let mut cliques_vec: Vec<Vec<u32>> = Vec::new();

            bron_kerbosch(&subgraph, &mut r, &mut p, &mut x, &mut cliques_vec);

            subclusters.extend(
                cliques_vec
                .into_iter()
                .map(|c| c.into_iter().collect()));
        }
    }

    let final_subclusters = connect_singletons_to_cliques(
        subclusters,
        ani_map,
        ani_cutoff,
        aligned_frac,
        id_to_name,
        bin_qualities,
        isolate_genomes,
        no_reassembly,
    );

    final_subclusters
}

// Detect maximal cliques
fn bron_kerbosch(
    graph: &Graph<u32, (), petgraph::Undirected>,
    r: &mut Vec<NodeIndex>,
    p: &mut Vec<NodeIndex>,
    x: &mut Vec<NodeIndex>,
    cliques: &mut Vec<Vec<u32>>,
) {
    if p.is_empty() && x.is_empty() {
        cliques.push(r.iter().map(|&n| graph[n]).collect());
        return;
    }

    // Tomita pivot: pick u in P∪X maximising |N(u) ∩ P| to minimise recursive calls
    let p_set: HashSet<NodeIndex> = p.iter().copied().collect();
    let pivot = p
        .iter()
        .chain(x.iter())
        .copied()
        .max_by_key(|&u| graph.neighbors(u).filter(|n| p_set.contains(n)).count())
        .unwrap();
    let pivot_neighbors: HashSet<NodeIndex> = graph.neighbors(pivot).collect();

    // Only recurse on vertices not adjacent to the pivot
    let candidates: Vec<NodeIndex> = p
        .iter()
        .copied()
        .filter(|v| !pivot_neighbors.contains(v))
        .collect();

    for v in candidates {
        r.push(v);

        // Build v's neighbor set once — O(1) membership checks replace O(degree) contains_edge
        let v_neighbors: HashSet<NodeIndex> = graph.neighbors(v).collect();
        let mut p_next: Vec<NodeIndex> = p.iter().copied().filter(|u| v_neighbors.contains(u)).collect();
        let mut x_next: Vec<NodeIndex> = x.iter().copied().filter(|u| v_neighbors.contains(u)).collect();

        bron_kerbosch(graph, r, &mut p_next, &mut x_next, cliques);
        r.pop();
        // move v from P to X
        if let Some(pos) = p.iter().position(|&u| u == v) {
            p.swap_remove(pos);
        }
        x.push(v);
    }
}

fn build_adj(
    component: &HashSet<u32>,
    ani_map: &HashMap<(u32, u32), AniData>,
    ani_cutoff: f32,
    aligned_frac: f32,
) -> HashMap<u32, HashSet<u32>> {
    let mut adj: HashMap<u32, HashSet<u32>> = HashMap::new();
    for &id in component {
        adj.entry(id).or_default();
    }

    let ids: Vec<u32> = component.iter().copied().collect();

    // Collect valid edges in parallel, then insert sequentially
    let edges: Vec<(u32, u32)> = (0..ids.len())
        .into_par_iter()
        .flat_map(|i| {
            let id1 = ids[i];
            ((i + 1)..ids.len())
                .filter_map(|j| {
                    let id2 = ids[j];
                    let key = if id1 <= id2 { (id1, id2) } else { (id2, id1) };
                    let data = ani_map.get(&key)?;
                    (data.ani >= ani_cutoff
                        && data.af_ref >= aligned_frac
                        && data.af_query >= aligned_frac)
                        .then_some((id1, id2))
                })
                .collect::<Vec<_>>()
        })
        .collect();

    for (id1, id2) in edges {
        adj.get_mut(&id1).unwrap().insert(id2);
        adj.get_mut(&id2).unwrap().insert(id1);
    }

    adj
}

fn extract_core_chunk(
    remaining: &mut HashSet<u32>,
    adj: &HashMap<u32, HashSet<u32>>,
    max_size: usize,
) -> HashSet<u32> {
    // Pick a seed: e.g. highest-degree node
    let &seed = remaining
        .iter()
        .max_by_key(|id| adj.get(id)
        .map_or(0, |neighs| neighs.len()))
        .expect("remaining non-empty");

    let mut core = HashSet::new();
    let mut frontier = vec![seed];
    core.insert(seed);

    while core.len() < max_size && 
        !frontier.is_empty() {

        let mut next_frontier = Vec::new();

        for node in frontier {
            if let Some(neighs) = adj.get(&node) {
                for &n in neighs {
                    if remaining.contains(&n) && !core.contains(&n) {
                        core.insert(n);
                        next_frontier.push(n);
                        if core.len() >= max_size {
                            break;
                        }
                    }
                }
            }
            if core.len() >= max_size {
                break;
            }
        }

        frontier = next_frontier;
    }

    // Remove core nodes from remaining
    for id in &core {
        remaining.remove(id);
    }

    core
}

fn build_subgraph_for_ids(
    core: &HashSet<u32>,
    adj: &HashMap<u32, HashSet<u32>>,
) -> Graph<u32, (), Undirected> {
    let mut subgraph: Graph<u32, (), Undirected> = Graph::default();
    let mut node_map: HashMap<u32, NodeIndex> = HashMap::new();

    for &id in core {
        let idx = subgraph.add_node(id);
        node_map.insert(id, idx);
    }

    for &id1 in core {
        if let Some(neighs) = adj.get(&id1) {
            for &id2 in neighs {
                if id1 < id2 && core.contains(&id2) {
                    let n1 = node_map[&id1];
                    let n2 = node_map[&id2];
                    subgraph.add_edge(n1, n2, ());
                }
            }
        }
    }
    subgraph
}

fn connect_singletons_to_cliques(
    clusters: Vec<HashSet<u32>>,
    ani_map: &HashMap<(u32, u32), AniData>,
    ani_cutoff: f32,
    aligned_frac: f32,
    id_to_name: &[String],
    bin_qualities: &HashMap<String, BinQuality>,
    isolate_genomes: &HashSet<String>,
    no_reassembly: bool,
) -> Vec<HashSet<u32>> {

    // Split into multi-node cliques & singleton nodes
    let mut cliques: Vec<HashSet<u32>> = Vec::new();
    let mut singletons: Vec<u32> = Vec::new();

    for cluster in clusters {
        if cluster.len() == 1 {
            if let Some(&only) = cluster.iter().next() {
                singletons.push(only);
            }
        } else {
            cliques.push(cluster);
        }
    }

    // Return (id, quality_score, is_isolate) for the best bin in a clique.
    // Isolate bins are always preferred over non-isolates; ties broken by quality score.
    let best_bin_of = |clique: &HashSet<u32>| -> (u32, f32, bool) {
        clique
            .iter()
            .copied()
            .map(|id| {
                let name = &id_to_name[id as usize];
                let score = bin_qualities.get(name).map_or(0.0, |q| q.score());
                let is_isolate = isolate_genomes.contains(name);
                (id, score, is_isolate)
            })
            .max_by(|(_, s1, i1), (_, s2, i2)| match (i1, i2) {
                (true, false) => std::cmp::Ordering::Greater,
                (false, true) => std::cmp::Ordering::Less,
                _ => s1.total_cmp(s2),
            })
            .unwrap_or((*clique.iter().next().unwrap(), 0.0, false))
    };

    // Returns true if (a_score, a_is_isolate) beats (b_score, b_is_isolate),
    // mirroring the priority rule used in MWIDS: isolate > non-isolate, then quality.
    let beats = |a_score: f32, a_iso: bool, b_score: f32, b_iso: bool| -> bool {
        (a_iso && !b_iso) || (a_iso == b_iso && a_score > b_score)
    };

    let mut leftover_singletons: Vec<HashSet<u32>> = Vec::new();

    for node in singletons {
        let node_name = &id_to_name[node as usize];
        let node_score = bin_qualities.get(node_name).map_or(0.0, |q| q.score());
        let node_is_isolate = isolate_genomes.contains(node_name);

        // A clique qualifies when its best representative has a passing ANI + AF
        // link to the query singleton.
        // Each entry: (clique_index, ani_to_best_bin, best_bin_score, best_bin_is_isolate)
        let mut qualified_cliques: Vec<(usize, f32, f32, bool)> = Vec::new();

        for (i, clique) in cliques.iter().enumerate() {
            let (best_id, best_score, best_is_isolate) = best_bin_of(clique);
            let key = if node <= best_id { (node, best_id) } else { (best_id, node) };

            let mut all_ok = true;

            for &member in clique.iter() {

                let key_m = if node <= member { (node, member) } else { (member, node) };

                let (ani, af_r, af_q) = ani_map.get(&key_m)
                    .map(|d| (d.ani, d.af_ref, d.af_query))
                    .unwrap_or((0.0, 0.0, 0.0));

                if ani < ani_cutoff || af_r < aligned_frac || af_q < aligned_frac {
                    all_ok = false;
                    break;
                }
            }

            if !all_ok {
                break;
            }

            if let Some(d) = ani_map.get(&key) {
                if d.ani >= ani_cutoff && d.af_ref >= aligned_frac && d.af_query >= aligned_frac {
                    qualified_cliques.push((i, d.ani, best_score, best_is_isolate));
                }
            }
        }

        match qualified_cliques.len() {
            0 => {
                // No clique representative links to this singleton — independent cluster.
                leftover_singletons.push(HashSet::from([node]));
            }
            1 => {
                // Exactly one potential clique — join it directly.
                cliques[qualified_cliques[0].0].insert(node);
            }
            _ => {
                // Multiple potential cliques. The query merges them all when it beats
                // every clique's best representative; otherwise it joins the closest one.
                if no_reassembly{
                    let query_is_best = qualified_cliques.iter().all(|&(_, _, best_score, best_iso)| {
                        beats(node_score, node_is_isolate, best_score, best_iso)
                    });

                    if query_is_best {
                        // Query is the highest-quality hub — merge all potential cliques
                        // into one and add the query.
                        let mut indices: Vec<usize> =
                            qualified_cliques.iter().map(|&(i, _, _, _)| i).collect();
                        // Process highest indices first so swap_remove doesn't
                        // invalidate the remaining indices.
                        indices.sort_unstable_by(|a, b| b.cmp(a));
                        let mut merged: HashSet<u32> = HashSet::from([node]);
                        for idx in indices {
                            merged.extend(cliques.swap_remove(idx));
                        }
                        cliques.push(merged);
                    } else {
                        // Query is not the best — assign it to the clique whose
                        // representative has the highest ANI to the query.
                        // Ties are broken by clique size (larger is better).
                        let best_idx = qualified_cliques
                            .iter()
                            .max_by(|&(ia, ani_a, _, _), &(ib, ani_b, _, _)| {
                                ani_a
                                    .total_cmp(ani_b)
                                    .then_with(|| cliques[*ia].len().cmp(&cliques[*ib].len()))
                            })
                            .map(|&(i, _, _, _)| i)
                            .unwrap();
                        cliques[best_idx].insert(node);
                    }
                } else {
                    // attaches to ALL qualifying cliques
                    for &(idx, _, _, _) in &qualified_cliques {
                        cliques[idx].insert(node);
                    }
                }
            }
        }
    }

    cliques.extend(leftover_singletons);
    cliques
}
