use crate::assess::BinQuality;
use crate::merge;
use crate::utility;
use log::{debug, info};
use petgraph::prelude::*;
use petgraph::visit::IntoNodeIdentifiers;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

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
            let key = if v_id <= u_id {
                (v_id, u_id)
            } else {
                (u_id, v_id)
            };
            // For an existing edge, ANI should exist; if not, treat as 0 contribution.
            if let Some(&ani) = ani_details.get(&key) {
                // edge exists iff ani >= c, so weight is non-negative
                let w = ani - ani_cutoff;
                if w > 0.0 {
                    sum += w;
                }
            }
        }
        wdeg.insert(v_id, sum);
    }
    wdeg
}

pub fn select_highconnectivity_bins(
    graph: &Graph<u32, (), Undirected>,
    ani_details: &HashMap<(u32, u32), f32>,
    ani_cutoff: f32,
    id_to_name: &Vec<String>,
    id_to_node: &HashMap<u32, NodeIndex>,
    bin_qualities: &HashMap<String, BinQuality>,
    result_dir: &Path,
) -> HashSet<String> {
    let isolate_genomes = HashSet::new();
    let (rep_set, rep_members) = select_highconnectivity_bins_with_memberships(
        graph,
        ani_details,
        ani_cutoff,
        id_to_name,
        id_to_node,
        bin_qualities,
        &isolate_genomes,
    );

    let _ = utility::write_membership_file(&rep_members, &result_dir.join("memberships.tsv"));
    let _ = write_quality_file(&rep_set, bin_qualities, result_dir);
    rep_set
}

pub fn select_highconnectivity_bins_with_memberships(
    graph: &Graph<u32, (), Undirected>,
    ani_details: &HashMap<(u32, u32), f32>,
    ani_cutoff: f32,
    id_to_name: &[String],
    id_to_node: &HashMap<u32, NodeIndex>,
    bin_qualities: &HashMap<String, BinQuality>,
    isolate_genomes: &HashSet<String>,
) -> (HashSet<String>, HashMap<String, Vec<String>>) {
    let connected_components = merge::compute_connected_components(&graph);

    let node_degrees: HashMap<u32, f32> = weighted_degree(&graph, &ani_details, ani_cutoff);

    let mut rep_set: HashSet<String> = HashSet::new();
    let mut assigned_rep: HashMap<u32, u32> = HashMap::new();

    for component in connected_components {
        if component.len() == 1 {
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
            let mut best_is_isolate = false;

            for &n in &component {
                if uncovered.contains(&n) && !blocked.contains(&n) {
                    // node degree includes edges from all nodes, including low quality bins
                    let wdeg = *node_degrees.get(&n).unwrap_or(&0.0);
                    let bin_name = &id_to_name[n as usize];
                    let bin_quality = &bin_qualities[bin_name];
                    let is_isolate = isolate_genomes.contains(bin_name);
                    let completeness = bin_quality.completeness;
                    let contamination = bin_quality.contamination;
                    let qscore = completeness - (5.0 * contamination);
                    debug!(
                        "Evaluating bin {}: wdeg {}, qscore {}",
                        bin_name, wdeg, qscore
                    );
                    let better = (is_isolate && !best_is_isolate)
                        || (is_isolate == best_is_isolate
                            && ((wdeg > best_wdeg)
                                || (wdeg == best_wdeg && qscore > best_quality)
                                || (wdeg == best_wdeg
                                    && qscore == best_quality
                                    && best_bin.map_or(true, |b| n < b))));

                    if better {
                        best_wdeg = wdeg;
                        best_quality = qscore;
                        best_is_isolate = is_isolate;
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

    let rep_members = rep_members_from_assigned(&assigned_rep, id_to_name);
    (rep_set, rep_members)
}

fn rep_members_from_assigned(
    assigned_rep: &HashMap<u32, u32>,
    id_to_name: &[String],
) -> HashMap<String, Vec<String>> {
    let mut rep_members: HashMap<String, Vec<String>> = HashMap::new();
    for (&member_id, &rep_id) in assigned_rep {
        let rep_name = id_to_name[rep_id as usize].clone();
        let member_name = id_to_name[member_id as usize].clone();
        rep_members.entry(rep_name.clone()).or_default();
        if member_id != rep_id {
            rep_members.entry(rep_name).or_default().push(member_name);
        }
    }
    rep_members
}

fn write_quality_file(
    rep_set: &HashSet<String>,
    bin_qualities: &HashMap<String, BinQuality>,
    result_dir: &Path,
) -> Result<(), Box<dyn std::error::Error>> {
    let output_file_path = result_dir.join("bins_checkm2_qualities.tsv");
    let output_file = File::create(&output_file_path)?;
    let mut writer = BufWriter::new(output_file);
    let mut reps: Vec<&String> = rep_set.iter().collect();
    reps.sort_unstable();

    writeln!(writer, "#Bin\tCompleteness\tContamination")?;
    for rep_name in reps {
        if let Some(quality) = bin_qualities.get(rep_name) {
            writeln!(
                writer,
                "{}\t{}\t{}",
                rep_name, quality.completeness, quality.contamination
            )?;
        }
    }
    info!(
        "Quality values of bins are written to {:?}",
        output_file_path
    );

    Ok(())
}