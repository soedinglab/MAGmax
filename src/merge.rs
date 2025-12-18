use bio::io::fasta;
use std::io::{BufWriter, Write};
use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};
use std::fs::{remove_file, File};
use std::io::{self, BufRead, BufReader};
use petgraph::graph::{Graph};
use petgraph::graph::NodeIndex;
use petgraph::visit::Dfs;
use petgraph::{Undirected};
use std::process::Command as ProcessCommand;
use log::{debug, error, info, warn};
use glob::glob;
use crate::assess::BinQuality;
use crate::cliques;

/// Compute all-vs-all ANI among bins
pub fn calc_ani(
    bins: &PathBuf,
    bin_qualities: &HashMap<String, BinQuality>,
    format: &String,
    anifile: Option<PathBuf>,
    ani_cutoff: f64,
    contamination_cutoff: f64,
    alignedfrac: f64,
    threads: usize
) -> Result<(Graph<u32, (), petgraph::Undirected>, HashMap<(u32, u32), f64>, Vec<String>), io::Error> {
    
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
                return Err(io::Error::new(io::ErrorKind::NotFound, "No fasta files found"));
            }

            let _ = get_ani(bin_files, &ani_output, threads)?;
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
            return Err(io::Error::new(io::ErrorKind::NotFound, "No fasta files found"));
        }
        
        let _ = get_ani(bin_files, &ani_output, threads)?;
    }
    // let mut bin_name_to_node: HashMap<String, prelude::NodeIndex> = HashMap::new();
    let mut bin_name_to_id: HashMap<String, u32> = HashMap::new();
    let mut id_to_name: Vec<String> = Vec::new();
    let mut id_to_node: HashMap<u32, NodeIndex> = HashMap::new();

    let file: File = File::open(ani_output.clone())?;
    let reader: BufReader<File> = io::BufReader::new(file);

    let mut graph: Graph<u32, (), Undirected> = Graph::default();

    // Add nodes to graph that pass quality filters
    for (bin, q) in bin_qualities {
        if q.contamination < contamination_cutoff && q.completeness > 20.0 {
            let id = get_or_assign_id(bin, &mut bin_name_to_id, &mut id_to_name);
            let node = graph.add_node(id);
            id_to_node.insert(id, node);
        }
    }

    let mut ani_details = HashMap::<(u32, u32), f64>::new();

    // Create a graph by add edge when ANI > ANI_threshold
    // When file is empty, no edge is formed and all nodes will be Singleton clusters.   
    for line in reader.lines().skip(1) {
        let line = line?;
        let columns: Vec<&str> = line.split('\t').collect();
        if columns.len() < 3 {
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

        let id1 = get_or_assign_id(
            &bin1, &mut bin_name_to_id, &mut id_to_name);
        let id2 = get_or_assign_id(
            &bin2, &mut bin_name_to_id, &mut id_to_name);
    
        let ani: f64 = columns[2].parse().unwrap_or(0.0);
        let alignfrac_ref: f64 = columns[3].parse().unwrap_or(0.0);
        let alignfrac_que: f64 = columns[4].parse().unwrap_or(0.0);
    
        // Skani reports pairs only if ANI is >= 80%
        if ani < ani_cutoff
            && alignfrac_ref >= alignedfrac
            && alignfrac_que >= alignedfrac {
            continue;
        }
    
        if let (Some(&node1), Some(&node2)) = 
            (id_to_node.get(&id1), id_to_node.get(&id2)) {
            graph.add_edge(node1, node2, ());
            let key = if id1 <= id2 { (id1, id2) } else { (id2, id1) };
            ani_details.insert(key, ani as f64);
        }
    }

    // Remove skani output file
    // remove_file(&ani_output).ok();
   Ok((graph, ani_details, id_to_name))
}

/// Find single-linkage connected components
pub fn get_connected_samples(
    graph: &Graph<u32, (), Undirected>,
    ani_details: &HashMap<(u32, u32), f64>,
    ani_cutoff: f64,
    id_to_name: &[String]
) -> Vec<HashSet<String>> {
    let mut visited = HashSet::new();
    let mut connected_components = Vec::new();

    for node_index in graph.node_indices() {
        if !visited.contains(&node_index) {
            // Start a new component
            let mut component = HashSet::new();
            let mut dfs = Dfs::new(&graph, node_index);

            while let Some(nx) = dfs.next(&graph) {
                if visited.insert(nx) {
                    let node_name = graph[nx].clone();
                    component.insert(node_name);
                }
            }
            connected_components.push(component);
        }
    }
    let mut connected_samples: Vec<HashSet<String>> = vec![];
    for component in connected_components {
        if component.len() <=2 {
            let component_names = component
                .into_iter()
                .map(|id| id_to_name[id as usize].clone())
                .collect();
            connected_samples.push(component_names);
        } else {
            let mut subclusters = 
                cliques::split_component_into_cliques(component, ani_details, ani_cutoff);
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
    format: &String,
) -> io::Result<()> {
    // Combine bins fasta into a single file
    let mut output_writer = File::create(
        combined_bins
        .join("combined.fasta"))?;
    for bin_samplename in bin_samplenames {
        let bin_file_path = inputdir.join(format!("{}.{}", bin_samplename, format));
        if bin_file_path.exists() {
            let bin_file = File::open(&bin_file_path)?;
            let reader = fasta::Reader::new(bin_file);

            for record in reader.records() {
                let record = record?;  // Get the record
                writeln!(output_writer, ">{}", format!("{}",record.id()))?;
                writeln!(output_writer, "{}", String::from_utf8_lossy(record.seq()))?;
            }
        } else {
            error!("Warning: File for bin '{}' does not exist at {:?}", bin_samplename, bin_file_path);
        }
    }
    Ok(())
}

/// Dereplicate final bins to remove any redundant bins
pub fn drep_finalbins(
    result_dir: &PathBuf,
    bin_qualities: &HashMap<String, BinQuality>,
    ani_cutoff: f64,
    threads: usize,
) -> Result<(), Box<dyn std::error::Error>> {
    let ani_output: PathBuf = result_dir.join("ani_edges");

    let finalbin_files: Vec<PathBuf> = glob(&format!("{}/*.fasta", result_dir.display()))
        .expect("Failed to read glob pattern")
        .filter_map(Result::ok)
        .collect();

    if let Err(e) = get_ani(
        finalbin_files.iter().map(|p| p.to_string_lossy().into_owned()).collect(), 
        &ani_output,
        threads,
    ) {
        error!("Error running skani for dereplication: {}", e);
        return Err(Box::new(e));
    };

    let file: File = File::open(ani_output.clone())?;
    let reader: BufReader<File> = io::BufReader::new(file);

    let mut bins_to_remove: HashSet<String> = HashSet::new();

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
    
        let ani: f64 = columns[2].parse().expect("Failed to parse ANI value as float from column 3");

        // Skani gives results for 80% aligned pairs
        if ani >= ani_cutoff {
            if let (Some(q1), Some(q2)) = (bin_qualities.get(&bin1), bin_qualities.get(&bin2)) {
                let score1 = q1.completeness - (5.0 * q1.contamination);
                let score2 = q2.completeness - (5.0 * q2.contamination);
    
                let worse_bin = if score1 > score2 {
                    &bin2
                } else if score1 < score2 {
                    &bin1
                } else if q1.contamination < q2.contamination {
                    &bin2
                } else {
                    &bin1
                };
                bins_to_remove.insert(worse_bin.clone());
            }
        }
    }

    debug!("Length of list with bins to remove: {}", bins_to_remove.len());
    let bin_names: HashSet<String> = finalbin_files
        .iter()
        .filter_map(|file| file.file_stem()?.to_str().map(|s| s.to_string()))
        .collect();

    let filtered_bin_names: HashSet<String> = bin_names
        .difference(&bins_to_remove)
        .cloned()
        .collect();

    // Remove redundant bins
    for bin in &bins_to_remove {
        let bin_file_path = result_dir.join(format!("{}.fasta", bin));
    
        if bin_file_path.exists() {
            remove_file(&bin_file_path).ok();
        }
    }
    
    if !cfg!(debug_assertions) {
        if let Err(e) = remove_file(&ani_output) {
            warn!("Failed to delete folder {:?}: {}", ani_output, e);
        }
    }
    
    // Write quality measures of bins
    let output_file_path = result_dir.join("bins_checkm2_qualities.tsv");
    let output_file = File::create(&output_file_path)?;
    let mut writer = BufWriter::new(output_file);

    writeln!(writer, "#Bin\tCompleteness\tContamination")?;
    let mut buffer = String::with_capacity(1024 * 1024);
    for (bin, quality) in bin_qualities.iter() {
        if filtered_bin_names.contains(bin) {
            buffer.push_str(&format!("{}\t{}\t{}\n", bin, quality.completeness, quality.contamination));
        }
    }
    writer.write_all(buffer.as_bytes())?;
    info!("Quality values of bins are written to {:?}", output_file_path);
    Ok(())
}

// Run skani
fn get_ani (
    inputbins:Vec<String>,
    ani_output: &PathBuf,
    threads: usize,
) -> Result<(), io::Error> {
    let mut command = ProcessCommand::new("skani");
    command.arg("triangle");
    command.args(&inputbins);
    command.arg("-E");
    command.arg("-t");
    command.arg(threads.to_string());
    
    if which::which("skani").is_err() {
        return Err(io::Error::new(io::ErrorKind::NotFound, "`skani` not found in PATH"));
    }
    
    let output = command.output()?;
    
    if !output.status.success() {
        return Err(io::Error::new(io::ErrorKind::Other, "skani triangle failed"));
    }
    std::fs::write(ani_output, output.stdout)?;
    Ok(())
}


fn get_or_assign_id(
    name: &str,
    map: &mut HashMap<String, u32>,
    names: &mut Vec<String>,
) -> u32 {
    if let Some(&id) = map.get(name) {
        return id;
    }
    let id = names.len() as u32;
    map.insert(name.to_owned(), id);
    names.push(name.to_owned());
    id
}