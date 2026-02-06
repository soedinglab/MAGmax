

use std::collections::{HashMap, HashSet};
use petgraph::prelude::*;
use petgraph::visit::{IntoNodeIdentifiers};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use crate::assess::BinQuality;
use crate::merge;
use log::{debug};

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
    result_dir: &Path,
    ) -> HashSet<String> {
    let connected_components = merge::compute_connected_components(&graph);

    let node_degrees: HashMap<u32, f32> = weighted_degree(&graph, &ani_details, ani_cutoff);
    
    let mut rep_set: HashSet<String> = HashSet::new();
    let mut assigned_rep: HashMap<u32, u32> = HashMap::new();

    for component in connected_components {
        
        if component.len()  == 1 {
            let bin_id = *component.iter().next().unwrap();
            let bin_name = &id_to_name[bin_id as usize];
            assigned_rep.insert(bin_id, bin_id);
            rep_set.insert(bin_name.clone());
            continue;
        }

        // nodes still needing domination
        let mut uncovered: HashSet<u32> = component.iter().copied().collect();
        // cannot be reps (adjacent to existing reps)
        let mut blocked: HashSet<u32> = HashSet::new();
        let mut chosen_reps: Vec<u32> = Vec::new();

        while !uncovered.is_empty() {
            // Select the best bin
            let mut best_bin: Option<u32> = None;
            let mut best_wdeg: f32 = -1.0;
            let mut best_quality: f32 = -1.0;

            for &n in &component {
                if uncovered.contains(&n) && !blocked.contains(&n) {
                    // node degree includes edges from all nodes, including low quality bins
                    let wdeg = *node_degrees.get(&n).unwrap_or(&0.0);
                    let bin_name = &id_to_name[n as usize];
                    let bin_quality = &bin_qualities[bin_name];
                    let completeness = bin_quality.completeness;
                    let contamination = bin_quality.contamination;
                    let qscore = completeness - (5.0 * contamination);
                    debug!("Evaluating bin {}: wdeg {}, qscore {}", bin_name, wdeg, qscore);
                    let better = (wdeg > best_wdeg)
                        || (wdeg == best_wdeg && qscore > best_quality)
                        || (wdeg == best_wdeg && qscore == best_quality
                        && best_bin.map_or(true, |b| n < b));

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

        // Representatives are chosen by max-weighted independent dominating set. 
        // Members are assigned to representative with highest connectivity, not by ANI.
        for rep_id in chosen_reps {
            rep_set.insert(id_to_name[rep_id as usize].clone());
        }

    }

    let _ = write_membershipfile(&assigned_rep, id_to_name, bin_qualities, result_dir);
    rep_set
}


fn write_membershipfile(
    assigned_rep: &HashMap<u32, u32>, // node_id -> rep_id
    id_to_name: &[String],
    bin_qualities: &HashMap<String, BinQuality>,
    result_dir: &Path
) -> Result<(), Box<dyn std::error::Error>> {   // rep_name -> [member_names]

    let membership_filepath: PathBuf = result_dir.join("memberships.tsv");
    let membership_file = File::create(&membership_filepath)?;
    let mut w = BufWriter::new(membership_file);

    let output_file_path = result_dir.join("bins_checkm2_qualities.tsv");
    let output_file = File::create(&output_file_path)?;
    let mut writer = BufWriter::new(output_file);

    let mut pairs: Vec<(u32, u32)> = assigned_rep
        .iter()
        .map(|(&m, &r)| (r, m))
        .collect();
    pairs.sort_unstable();

    writeln!(writer, "bin_name\tcompleteness\tcontamination")?;
    let unique_reps: HashSet<u32> = assigned_rep.values().copied().collect();
    for &rep_id in unique_reps.iter() {
        let rep_name = &id_to_name[rep_id as usize];
        if let Some(quality) = bin_qualities.get(rep_name) {
            writeln!(
                writer,
                "{}\t{:.2}\t{:.2}",
                rep_name,
                quality.completeness,
                quality.contamination
            )?;
        }
    }

    let mut cur_rep: Option<u32> = None;
    let mut members: Vec<String> = Vec::new();
    writeln!(w, "#representative\tmember_genomes")?;

    for (rep_id, member_id) in pairs {
        if cur_rep != Some(rep_id) {
            // flush previous rep
            if let Some(r) = cur_rep {
                writeln!(
                    w,
                    "{}\t{}",
                    id_to_name[r as usize],
                    members.join(",")
                )?;
                members.clear();
            }
            cur_rep = Some(rep_id);
        }
        if member_id != rep_id {
            members.push(id_to_name[member_id as usize].clone());
        }
    }

    // flush last rep
    if let Some(r) = cur_rep {
        writeln!(
            w,
            "{}\t{}",
            id_to_name[r as usize],
            members.join(",")
        )?;
    }
    Ok(())
}