

use std::collections::{HashMap, HashSet};
use petgraph::prelude::*;
use petgraph::visit::{IntoNodeIdentifiers};
use crate::assess::BinQuality;
use crate::merge;



/// Greedy Max-Weighted Independent Dominating Set (MWIDS) functions

// Compute weighted degree for each node in the graph
fn weighted_degree(
    graph: &Graph<u32, (), Undirected>,
    ani_details: &HashMap<(u32, u32), f32>,
    ani_cutoff: f32,
) -> HashMap<u32, f32> {
    let mut wdeg: HashMap<u32, f32> = HashMap::new();

    for ni in graph.node_identifiers() {
        let v_id = graph[ni];
        let mut sum = 0.0f32;

        for nb in graph.neighbors(ni) {
            let u_id = graph[nb];
            let key = if v_id <= u_id { (v_id, u_id) } else { (u_id, v_id) };
            // For an existing edge, ANI should exist; if not, treat as 0 contribution.
            if let Some(&ani) = ani_details.get(&key) {
                // edge exists iff ani >= c, so weight is non-negative
                let w = ani - ani_cutoff;
                if w > 0.0 { sum += w; }
            }
        }
        wdeg.insert(v_id, sum);
    }
    wdeg
}

pub fn select_highconnectivity_bins(
    graph  : &Graph<u32, (), Undirected>,
    ani_details: &HashMap<(u32, u32), f32>,
    ani_cutoff: f32,
    id_to_name: &Vec<String>,
    id_to_node: &HashMap<u32, NodeIndex>,
    bin_qualities: &HashMap<String, BinQuality>,
    completeness_cutoff: f32,
    contamination_cutoff: f32,
    ) -> HashSet<String> {
    let connected_components = merge::compute_connected_components(&graph);

    let node_degrees: HashMap<u32, f32> = weighted_degree(&graph, &ani_details, ani_cutoff);
    
    let mut rep_set: HashSet<String> = HashSet::new();
    let mut assigned_rep: HashMap<u32, u32> = HashMap::new();

    for component in connected_components {
        let mut candidate_bins: Vec<u32> = Vec::new();
        for &node_id in &component {
            let bin_name = &id_to_name[node_id as usize];
            let bin_quality = &bin_qualities[bin_name];

            if bin_quality.completeness >= completeness_cutoff &&
               bin_quality.contamination <= contamination_cutoff {
                candidate_bins.push(node_id);
            }
        }

        if candidate_bins.len()  == 1 {
            let bin_id = candidate_bins[0];
            let bin_name = &id_to_name[bin_id as usize];
            println!("Selected bin: {}", bin_name);
            rep_set.insert(bin_name.clone());
            continue;
        }

        let mut uncovered: HashSet<u32> = component.clone(); // nodes still needing domination
        let mut blocked: HashSet<u32> = HashSet::new(); // cannot be reps (adjacent to existing reps)
        let mut chosen_reps: Vec<u32> = Vec::new();

        while !uncovered.is_empty() {
            // Select the best candidate
            let mut best_bin: Option<u32> = None;
            let mut best_wdeg: f32 = -1.0;
            let mut best_quality: f32 = -1.0;

            for &n in &candidate_bins {
                if uncovered.contains(&n) && !blocked.contains(&n) {
                    let wdeg = *node_degrees.get(&n).unwrap_or(&0.0);
                    let bin_name = &id_to_name[n as usize];
                    let bin_quality = &bin_qualities[bin_name];
                    let completeness = bin_quality.completeness;
                    let contamination = bin_quality.contamination;
                    let qscore = completeness - (5.0 * contamination);
                    let better = (wdeg > best_wdeg)
                        || (wdeg == best_wdeg && qscore > best_quality)
                        || (wdeg == best_wdeg && qscore == best_quality && best_bin.map_or(true, |b| n < b));

                    if better {
                        best_wdeg = wdeg;
                        best_quality = qscore;
                        best_bin = Some(n);
                    }
                }
            }
            
            let r_id = best_bin.unwrap_or_else(|| *uncovered.iter().next().unwrap());
            chosen_reps.push(r_id);

            let r_ni = id_to_node[&r_id];

            blocked.insert(r_id);
            for nb in graph.neighbors(r_ni) {
                let nb_id = graph[nb];
                if component.contains(&nb_id) {
                    blocked.insert(nb_id);
                }
            }

            if uncovered.remove(&r_id) {
                assigned_rep.entry(r_id).or_insert(r_id);
            }
            for nb in graph.neighbors(r_ni) {
                let nb_id = graph[nb];
                if component.contains(&nb_id) && uncovered.remove(&nb_id) {
                    assigned_rep.entry(nb_id).or_insert(r_id);
                }
            }
        }

         for rep_id in chosen_reps {
            rep_set.insert(id_to_name[rep_id as usize].clone());
        }

    }

    let rep_to_members_names = reps_to_members_by_name(&assigned_rep, &id_to_name);
    for (rep_name, member_names) in rep_to_members_names.iter() {
        println!("Representative bin: {}", rep_name);
        println!("  Member bins:");
        for member in member_names {
            println!("    {}", member);
        }
    }
    rep_set
}


fn reps_to_members_by_name(
    assigned_rep: &HashMap<u32, u32>, // node_id -> rep_id
    id_to_name: &[String],
) -> HashMap<String, Vec<String>> {   // rep_name -> [member_names]
    let mut rep_to_members: HashMap<String, Vec<String>> = HashMap::new();

    for (&node_id, &rep_id) in assigned_rep.iter() {
        let rep_name = id_to_name[rep_id as usize].clone();
        let member_name = id_to_name[node_id as usize].clone();
        rep_to_members.entry(rep_name).or_default().push(member_name);
    }

    // Optional: deterministic ordering
    for members in rep_to_members.values_mut() {
        members.sort();
    }

    rep_to_members
}
