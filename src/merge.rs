use crate::assess::BinQuality;
use crate::cliques;
use crate::utility;
use bio::io::fasta;
use glob::glob;
use log::{error, info, warn};
use petgraph::graph::Graph;
use petgraph::graph::NodeIndex;
use petgraph::visit::Dfs;
use petgraph::Undirected;
use std::collections::{HashMap, HashSet};
use std::fs::{remove_file, File};
use std::io::{self, BufRead, BufReader};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

/// Compute all-vs-all ANI among bins
pub fn calc_ani(
    bins: &Path,
    bin_qualities: &HashMap<String, BinQuality>,
    format: &str,
    anifile: Option<PathBuf>,
    ani_cutoff: f32,
    completeness_cutoff: f32,
    contamination_cutoff: f32,
    alignedfrac: f32,
    no_reassembly: bool,
    threads: usize,
) -> Result<
    (
        Graph<u32, (), petgraph::Undirected>,
        HashMap<(u32, u32), f32>,
        Vec<String>,
        HashMap<u32, NodeIndex>,
        HashMap<(u32, u32), f32>,
        HashMap<(u32, u32), f32>,
    ),
    io::Error,
> {
    let ani_output: PathBuf;

    if let Some(provided_path) = anifile {
        if provided_path.exists() {
            info!("Using provided ANI file at {:?}", provided_path);
            ani_output = provided_path;
        } else {
            info!("Provided ANI file not found; computing ANI ...");
            ani_output = bins.join("ani_edges");

            let bin_files: Vec<String> = glob(&format!("{}/*.{}", bins.display(), format))
                .expect("Failed to read glob pattern")
                .filter_map(Result::ok)
                .map(|path| path.to_string_lossy().into_owned())
                .collect();

            if bin_files.is_empty() {
                error!("No fasta files found in {:?}", bins);
                return Err(io::Error::new(
                    io::ErrorKind::NotFound,
                    "No fasta files found",
                ));
            }
            info!("Calculating ANI between bins using skani ...");
            get_ani(bin_files, &ani_output, threads)?;
        }
    } else {
        ani_output = bins.join("ani_edges");
        let bin_files: Vec<String> = glob(&format!("{}/*.{}", bins.display(), format))
            .expect("Failed to read glob pattern")
            .filter_map(Result::ok)
            .map(|path| path.to_string_lossy().into_owned())
            .collect();

        if bin_files.is_empty() {
            error!("No fasta files found in {:?}", bins);
            return Err(io::Error::new(
                io::ErrorKind::NotFound,
                "No fasta files found",
            ));
        }
        info!("Calculating ANI between bins using skani ...");
        get_ani(bin_files, &ani_output, threads)?;
    }

    let mut bin_name_to_id: HashMap<String, u32> = HashMap::new();
    let mut id_to_name: Vec<String> = Vec::new();
    let mut id_to_node: HashMap<u32, NodeIndex> = HashMap::new();

    let file: File = File::open(ani_output.clone())?;
    let reader: BufReader<File> = io::BufReader::new(file);

    let mut graph: Graph<u32, (), Undirected> = Graph::default();

    // Add nodes to graph that pass quality
    if no_reassembly {
        for (bin, q) in bin_qualities {
            // Filter bins based on both completeness and contamination cutoffs before constructing the graph
            if q.contamination < contamination_cutoff && q.completeness >= completeness_cutoff {
                let id = get_or_assign_id(bin, &mut bin_name_to_id, &mut id_to_name);
                let node = graph.add_node(id);
                id_to_node.insert(id, node);
            }
        }
    } else {
        for (bin, q) in bin_qualities {
            // In reassembly mode, filter contaminated bins only before constructing the graph
            if q.contamination < contamination_cutoff && q.completeness > 20.0 {
                let id = get_or_assign_id(bin, &mut bin_name_to_id, &mut id_to_name);
                let node = graph.add_node(id);
                id_to_node.insert(id, node);
            }
        }
    }

    let mut ani_details = HashMap::<(u32, u32), f32>::new();
    let mut af_ref = HashMap::<(u32, u32), f32>::new();
    let mut af_query = HashMap::<(u32, u32), f32>::new();

    // Create a graph by add edge when ANI > ANI_threshold
    // When file is empty, no edge is formed and all nodes will be Singleton clusters.
    for line in reader.lines().skip(1) {
        let line = line?;
        let columns: Vec<&str> = line.split('\t').collect();

        if columns.len() < 5 {
            continue;
        }

        let bin1 = Path::new(columns[0])
            .file_stem()
            .map(|name| name.to_string_lossy().into_owned())
            .unwrap_or_else(|| columns[0].to_string());

        let bin2 = Path::new(columns[1])
            .file_stem()
            .map(|name| name.to_string_lossy().into_owned())
            .unwrap_or_else(|| columns[1].to_string());

        let (Some(&id1), Some(&id2)) = (bin_name_to_id.get(&bin1), bin_name_to_id.get(&bin2))
        else {
            continue;
        };

        let ani: f32 = columns[2].parse().unwrap_or(0.0);
        let alignfrac_ref: f32 = columns[3].parse().unwrap_or(0.0);
        let alignfrac_que: f32 = columns[4].parse().unwrap_or(0.0);

        let key = if id1 <= id2 { (id1, id2) } else { (id2, id1) };
        ani_details.insert(key, ani as f32);
        af_ref.insert(key, alignfrac_ref as f32);
        af_query.insert(key, alignfrac_que as f32);

        // Skani reports pairs only if ANI is >= 80%
        if ani < ani_cutoff || alignfrac_ref < alignedfrac || alignfrac_que < alignedfrac {
            continue;
        }

        if let (Some(&node1), Some(&node2)) = (id_to_node.get(&id1), id_to_node.get(&id2)) {
            graph.add_edge(node1, node2, ());
        }
    }

    // Remove skani output file
    // remove_file(&ani_output).ok();
    Ok((graph, ani_details, id_to_name, id_to_node, af_ref, af_query))
}

/// pub fn single-linkage connected components

pub fn compute_connected_components(graph: &Graph<u32, (), Undirected>) -> Vec<HashSet<u32>> {
    let mut visited = HashSet::new();
    let mut connected_components = Vec::new();

    for node_index in graph.node_indices() {
        if !visited.contains(&node_index) {
            // Start a new component
            let mut component = HashSet::new();
            let mut dfs = Dfs::new(&graph, node_index);

            while let Some(nx) = dfs.next(&graph) {
                if visited.insert(nx) {
                    let node_id = graph[nx];

                    component.insert(node_id);
                }
            }
            connected_components.push(component);
        }
    }

    connected_components
}
/// Find single-linkage connected components
pub fn get_connected_samples(
    graph: &Graph<u32, (), Undirected>,
    ani_details: &HashMap<(u32, u32), f32>,
    ani_cutoff: f32,
    id_to_name: &[String],
    alignedfrac: f32,
    af_ref: &HashMap<(u32, u32), f32>,
    af_query: &HashMap<(u32, u32), f32>,
) -> Vec<HashSet<String>> {
    let connected_components = compute_connected_components(graph);

    let mut connected_samples: Vec<HashSet<String>> = Vec::new();
    for component in connected_components {
        if component.len() <= 2 {
            let component_names: HashSet<String> = component
                .into_iter()
                .map(|id| id_to_name[id as usize].clone())
                .collect();
            connected_samples.push(component_names);
        } else {
            let mut subclusters = cliques::split_component_into_cliques(
                component,
                ani_details,
                ani_cutoff,
                alignedfrac,
                af_ref,
                af_query,
            );
            for cluster in subclusters.drain(..) {
                let component_names: HashSet<String> = cluster
                    .into_iter()
                    .map(|id| id_to_name[id as usize].clone())
                    .collect();
                connected_samples.push(component_names);
            }
        }
    }
    connected_samples
}

/// Merged bin files within the cluster
pub fn combine_fastabins(
    inputdir: &Path,
    bin_samplenames: &HashSet<String>,
    combined_bins: &Path,
    format: &str,
) -> io::Result<()> {
    // Combine bins fasta into a single file

    let out = File::create(combined_bins.join("combined.fasta"))?;
    let mut output_writer = BufWriter::new(out);

    for bin_samplename in bin_samplenames {
        let bin_file_path = inputdir.join(format!("{}.{}", bin_samplename, format));
        if !bin_file_path.exists() {
            warn!(
                "Warning: File for bin '{}' does not exist at {:?}",
                bin_samplename, bin_file_path
            );
            continue;
        }

        let bin_file = File::open(&bin_file_path)?;
        let reader = fasta::Reader::new(bin_file);

        for record in reader.records() {
            let record = record?; // Get the record
            writeln!(output_writer, ">{}", record.id())?;
            writeln!(output_writer, "{}", String::from_utf8_lossy(record.seq()))?;
        }
    }
    Ok(())
}

/// Dereplicate final bins to remove any redundant bins
pub fn drep_finalbins(
    result_dir: &Path,
    bin_qualities: &HashMap<String, BinQuality>,
    ani_details: &HashMap<(u32, u32), f32>,
    id_to_name: &[String],
    af_ref: &HashMap<(u32, u32), f32>,
    af_query: &HashMap<(u32, u32), f32>,
    ani_cutoff: f32,
    alignedfrac: f32,
    threads: usize,
    noreassembly: bool,
    memberships_map: &HashMap<String, String>,
    format: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let finalbin_files: Vec<PathBuf> = glob(&format!("{}/*.{}", result_dir.display(), format))
        .expect("Failed to read glob pattern")
        .filter_map(Result::ok)
        .collect();

    let bin_names: HashSet<String> = finalbin_files
        .iter()
        .filter_map(|file| file.file_stem().map(|s| s.to_string_lossy().into_owned()))
        .collect();

    let mut bins_to_remove: HashSet<String> = HashSet::new();
    let mut bins_pair: HashMap<String, String> = HashMap::new();

    if noreassembly {
        for (&pair @ (id1, id2), &ani) in ani_details.iter() {
            // IDs -> names
            let Some(bin1) = id_to_name.get(id1 as usize) else {
                continue;
            };
            let Some(bin2) = id_to_name.get(id2 as usize) else {
                continue;
            };

            // Only consider pairs where both bins exist in final bins
            if !(bin_names.contains(bin1) && bin_names.contains(bin2)) {
                continue;
            }

            if ani < ani_cutoff {
                continue;
            }

            let af_r = af_ref
                .get(&pair)
                .or_else(|| af_ref.get(&(id2, id1)))
                .copied()
                .unwrap_or(0.0);

            let af_q = af_query
                .get(&pair)
                .or_else(|| af_query.get(&(id2, id1)))
                .copied()
                .unwrap_or(0.0);

            if af_r < alignedfrac || af_q < alignedfrac {
                continue;
            }

            let Some(q1) = bin_qualities.get(bin1) else {
                continue;
            };
            let Some(q2) = bin_qualities.get(bin2) else {
                continue;
            };

            let worse_bin = find_worsebin(bin1.as_str(), bin2.as_str(), q1, q2);
            let best_bin = if worse_bin == bin1 {
                bin2.as_str()
            } else {
                bin1.as_str()
            };

            bins_to_remove.insert(worse_bin.to_string());
            add_edge_keep_best(&mut bins_pair, worse_bin, best_bin, bin_qualities);
        }
    } else {
        let ani_output: PathBuf = result_dir.join("ani_edges");

        if let Err(e) = get_ani(
            finalbin_files
                .iter()
                .map(|p| p.to_string_lossy().into_owned())
                .collect(),
            &ani_output,
            threads,
        ) {
            error!("Error running skani for dereplication: {}", e);
            return Err(Box::new(e));
        };

        let file: File = File::open(ani_output.clone())?;
        let reader: BufReader<File> = io::BufReader::new(file);

        for line in reader.lines().skip(1) {
            let line: String = line?;
            let columns: Vec<&str> = line.split('\t').collect();

            let bin1 = Path::new(columns[0])
                .file_stem()
                .map(|name| name.to_string_lossy().into_owned())
                .unwrap_or_else(|| columns[0].to_string());

            let bin2 = Path::new(columns[1])
                .file_stem()
                .map(|name| name.to_string_lossy().into_owned())
                .unwrap_or_else(|| columns[1].to_string());

            let ani: f32 = columns[2]
                .parse()
                .expect("Failed to parse ANI value as float from column 3");
            let alignfrac_ref: f32 = columns[3].parse().unwrap_or(0.0);
            let alignfrac_que: f32 = columns[4].parse().unwrap_or(0.0);
            // Skani gives results for 80% aligned pairs
            if ani >= ani_cutoff && alignfrac_ref >= alignedfrac && alignfrac_que >= alignedfrac {
                if let (Some(q1), Some(q2)) = (bin_qualities.get(&bin1), bin_qualities.get(&bin2)) {
                    let worse_bin = find_worsebin(bin1.as_str(), bin2.as_str(), q1, q2);
                    let best_bin = if worse_bin == bin1 {
                        bin2.as_str()
                    } else {
                        bin1.as_str()
                    };
                    bins_to_remove.insert(worse_bin.to_string());
                    add_edge_keep_best(&mut bins_pair, worse_bin, best_bin, bin_qualities);
                }
            }
        }
        if !cfg!(debug_assertions) {
            if let Err(e) = remove_file(&ani_output) {
                warn!("Failed to delete file {:?}: {}", ani_output, e);
            }
        }
    }

    let updated_memberships = update_memberships_map(memberships_map, &bins_pair, &bins_to_remove);

    let filtered_bin_names: HashSet<String> =
        bin_names.difference(&bins_to_remove).cloned().collect();

    // Remove redundant bins
    for bin in &bins_to_remove {
        let bin_file_path = result_dir.join(format!("{}.{}", bin, &format));

        if bin_file_path.exists() {
            remove_file(&bin_file_path).ok();
        }
    }

    let membership_filepath: PathBuf = result_dir.join("memberships.tsv");
    let rep_members = utility::rep_members_from_member_rep(&updated_memberships);
    utility::write_membership_file(&rep_members, &membership_filepath)?;

    // Write quality measures of bins
    let output_file_path = result_dir.join("bins_checkm2_qualities.tsv");
    let output_file = File::create(&output_file_path)?;
    let mut writer = BufWriter::new(output_file);

    writeln!(writer, "#Bin\tCompleteness\tContamination")?;
    let mut buffer = String::with_capacity(1024 * 1024);
    for (bin, quality) in bin_qualities.iter() {
        if filtered_bin_names.contains(bin) {
            buffer.push_str(&format!(
                "{}\t{}\t{}\n",
                bin, quality.completeness, quality.contamination
            ));
        }
    }
    writer.write_all(buffer.as_bytes())?;
    info!(
        "Quality values of bins are written to {:?}",
        output_file_path
    );
    Ok(())
}

// Run skani
pub fn get_ani(inputbins: Vec<String>, ani_output: &PathBuf, threads: usize) -> Result<(), io::Error> {
    if which::which("skani").is_err() {
        return Err(io::Error::new(
            io::ErrorKind::NotFound,
            "`skani` not found in PATH",
        ));
    }

    let output_file = File::create(ani_output)?;
    let status = Command::new("skani")
        .arg("triangle")
        .args(&inputbins)
        .arg("-E")
        .arg("-t")
        .arg(threads.to_string())
        .stdout(Stdio::from(output_file))
        .stderr(Stdio::null())
        .status()?;

    if !status.success() {
        return Err(io::Error::new(
            io::ErrorKind::Other,
            "skani triangle failed",
        ));
    }

    Ok(())
}

fn get_or_assign_id(name: &str, map: &mut HashMap<String, u32>, names: &mut Vec<String>) -> u32 {
    if let Some(&id) = map.get(name) {
        return id;
    }
    let id = u32::try_from(names.len()).expect("too many names");
    map.insert(name.to_owned(), id);
    names.push(name.to_owned());
    id
}

// Get poorer quality bin
fn find_worsebin<'a>(bin1: &'a str, bin2: &'a str, q1: &BinQuality, q2: &BinQuality) -> &'a str {
    let score1 = q1.completeness - (5.0 * q1.contamination);
    let score2 = q2.completeness - (5.0 * q2.contamination);

    if score1 > score2 {
        bin2
    } else if score1 < score2 {
        bin1
    } else if q1.contamination < q2.contamination {
        // keep lower contamination, remove higher contamination
        bin2
    } else {
        bin1
    }
}

fn add_edge_keep_best(
    bins_pair: &mut HashMap<String, String>,
    worse: &str,
    better: &str,
    bin_qualities: &HashMap<String, BinQuality>,
) {
    match bins_pair.get(worse) {
        None => {
            bins_pair.insert(worse.to_string(), better.to_string());
        }
        Some(prev_better) if prev_better == better => {}
        Some(prev_better) => {
            let (Some(q_prev), Some(q_new)) =
                (bin_qualities.get(prev_better), bin_qualities.get(better))
            else {
                return;
            };

            let worse_of_candidates = find_worsebin(prev_better.as_str(), better, q_prev, q_new);
            let best_of_candidates = if worse_of_candidates == prev_better.as_str() {
                better
            } else {
                prev_better.as_str()
            };

            bins_pair.insert(worse.to_string(), best_of_candidates.to_string());
        }
    }
}

fn update_memberships_map(
    memberships_map: &HashMap<String, String>,
    bins_pair: &HashMap<String, String>,
    bins_to_remove: &HashSet<String>,
) -> HashMap<String, String> {
    // Memoization cache: rep -> final_rep
    let mut memo: HashMap<String, String> = HashMap::new();

    fn resolve_rep(
        rep: &str,
        bins_to_remove: &HashSet<String>,
        next: &HashMap<String, String>,
        memo: &mut HashMap<String, String>,
    ) -> Option<String> {
        if let Some(v) = memo.get(rep) {
            return Some(v.clone());
        }

        let mut cur = rep;
        let mut path: Vec<String> = Vec::new();

        while bins_to_remove.contains(cur) {
            path.push(cur.to_string());
            let Some(n) = next.get(cur) else {
                return None;
            };
            cur = n.as_str();
        }

        let final_rep = cur.to_string();

        for p in path {
            memo.insert(p, final_rep.clone());
        }
        memo.insert(rep.to_string(), final_rep.clone());

        Some(final_rep)
    }

    let mut updated = HashMap::with_capacity(memberships_map.len());
    for (member, rep) in memberships_map {
        let new_rep = resolve_rep(rep.as_str(), bins_to_remove, &bins_pair, &mut memo)
            .unwrap_or_else(|| rep.clone());

        updated.insert(member.clone(), new_rep);
    }

    updated
}
