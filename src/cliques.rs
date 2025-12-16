use std::collections::{HashMap, HashSet};
use petgraph::graph::{Graph};
use petgraph::graph::NodeIndex;
use petgraph::{Undirected};
use log::{debug};
use rayon::prelude::*;
use std::thread;

// Get clique clusters
pub fn split_component_into_cliques(
    component: HashSet<u32>,
    ani_details: &HashMap<(u32, u32), f64>,
    ani_cutoff: f64,
) -> Vec<HashSet<u32>> {

    let adj = build_adj(&component, ani_details, ani_cutoff);

    let mut remaining = component.clone();
    let mut subclusters: Vec<HashSet<u32>> = Vec::new();
    let mut cores: Vec<HashSet<u32>> = Vec::new();

    const MAX_CLIQUE_CORE_SIZE: usize = 80;
    while remaining.len() > MAX_CLIQUE_CORE_SIZE {
        if remaining.len() > 500{
            debug!("Remaining nodes to process for cliques: {:?}", remaining);
        }
        let core = 
            extract_core_chunk(&mut remaining, &adj, MAX_CLIQUE_CORE_SIZE);
        if core.len() <= 1 {
            // treat as singleton(s)
            for id in core {
                let mut s = HashSet::new();
                s.insert(id);
                subclusters.push(s);
            }
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
        // let subgraph = build_subgraph_for_ids(&core, &adj);
        // let mut r: Vec<NodeIndex> = Vec::new();
        // let mut p: Vec<NodeIndex> = subgraph.node_indices().collect();
        // let mut x: Vec<NodeIndex> = Vec::new();
        // let mut cliques_vec: Vec<Vec<u32>> = Vec::new();

        // bron_kerbosch(&subgraph, &mut r, &mut p, &mut x, &mut cliques_vec);

        // for c in cliques_vec {
        //     subclusters.push(c.into_iter().collect());
        // }
    }

    if !remaining.is_empty() {
        if remaining.len() == 1 {
            let mut s = HashSet::new();
            s.insert(*remaining.iter().next().unwrap());
            subclusters.push(s);
        } else {
            let subgraph = build_subgraph_for_ids(&remaining, &adj);

            let mut r: Vec<NodeIndex> = Vec::new();
            let mut p: Vec<NodeIndex> = subgraph.node_indices().collect();
            let mut x: Vec<NodeIndex> = Vec::new();
            let mut cliques_vec: Vec<Vec<u32>> = Vec::new();

            bron_kerbosch(&subgraph, &mut r, &mut p, &mut x, &mut cliques_vec);

            for c in cliques_vec {
                subclusters.push(c.into_iter().collect());
            }
        }
    }

    let subclusters = connect_singletons_to_cliques(subclusters.clone(), ani_details, ani_cutoff);
    
    subclusters
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

    let pivot = p.iter()
        .chain(x.iter()).next().copied().unwrap(); // Choose a pivot
    let neighbors: HashSet<NodeIndex> = graph
        .neighbors(pivot)
        .collect();

    let mut candidates: Vec<NodeIndex> = Vec::new();
    
    for &v in p.iter() {
        if !neighbors.contains(&v) {
            candidates.push(v);
        }
    }
    
    for v in candidates {
        r.push(v);

        let mut p_next: Vec<NodeIndex> = Vec::new();

        for &u in p.iter() {
            if graph.contains_edge(v, u) {
                p_next.push(u);
            }
        }
        let mut x_next: Vec<NodeIndex> = Vec::new();
        for &u in x.iter() {
            if graph.contains_edge(v, u) {
                x_next.push(u);
            }
        }

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
    ani_details: &HashMap<(u32, u32), f64>,
    ani_cutoff: f64,
) -> HashMap<u32, HashSet<u32>> {
    let mut adj: HashMap<u32, HashSet<u32>> = HashMap::new();

    for &id in component {
        adj.entry(id).or_default();
    }

    let ids: Vec<u32> = component.iter().copied().collect();
    for i in 0..ids.len() {
        for j in (i + 1)..ids.len() {
            let id1 = ids[i];
            let id2 = ids[j];

            let key = if id1 <= id2 { (id1, id2) } else { (id2, id1) };

            if let Some(&ani) = ani_details.get(&key) {
                if ani >= ani_cutoff {
                    adj.get_mut(&id1).unwrap().insert(id2);
                    adj.get_mut(&id2).unwrap().insert(id1);
                }
            }
        }
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
        .max_by_key(|id| adj.get(id).map_or(0, |neighs| neighs.len()))
        .expect("remaining non-empty");

    let mut core = HashSet::new();
    let mut frontier = vec![seed];
    core.insert(seed);

    while core.len() < max_size {
        if frontier.is_empty() {
            break;
        }

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
    debug!("Building subgraph for core of size {} using thread {:?}", core.len(), thread::current().id());
    for &id1 in core {
        if let Some(neighs) = adj.get(&id1) {
            for &id2 in neighs {
                if id1 < id2 && core.contains(&id2) {
                    if let (Some(&n1), Some(&n2)) = (node_map.get(&id1), node_map.get(&id2)) {
                        subgraph.add_edge(n1, n2, ());
                    }
                }
            }
        }
    }

    subgraph
}

fn connect_singletons_to_cliques(
    clusters: Vec<HashSet<u32>>,
    ani_details: &HashMap<(u32, u32), f64>,
    ani_cutoff: f64,
) -> Vec<HashSet<u32>> {

    // Split into multi-node cliques & singleton nodes
    let mut cliques: Vec<HashSet<u32>> = Vec::new();
    let mut singletons: Vec<u32> = Vec::new();

    for cluster in clusters {
        if cluster.len() == 1 {
            singletons.push(*cluster.iter().next().unwrap());
        } else {
            cliques.push(cluster);
        }
    }

    let mut leftover_singletons: Vec<HashSet<u32>> = Vec::new();

    for node in singletons {
        let mut qualified_cliques = Vec::new();

        // Check which cliques this node can join
        for (i, clique) in cliques.iter().enumerate() {
            let mut all_ok = true;

            for &member in clique {
                let (a, b) = if node <= member { (node, member) } else { (member, node) };

                match ani_details.get(&(a, b)) {
                    Some(&ani) if ani >= ani_cutoff => { /* ok */ }
                    _ => {
                        all_ok = false;
                        break;
                    }
                }
            }

            if all_ok {
                qualified_cliques.push(i);
            }
        }

        match qualified_cliques.len() {
            0 => {
                // stays singleton
                leftover_singletons.push(HashSet::from([node]));
            }
            1 => {
                // exactly one match → glue into that one clique
                let idx = qualified_cliques[0];
                cliques[idx].insert(node);
            }
            _ => {
                // attaches to ALL qualifying cliques
                for &idx in &qualified_cliques {
                    cliques[idx].insert(node);
                }
            }
        }
    }

    // return cliques + leftover singletons
    cliques.extend(leftover_singletons);
    cliques
}
