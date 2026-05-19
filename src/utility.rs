use log::{error, info};
use rayon::prelude::*;
use rayon::ThreadPool;
use std::collections::{HashMap, HashSet};
use std::fs::{self, read_to_string, File};
use std::io::{self, BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::exit;
use std::sync::{Arc, Mutex};

// Helper function
pub fn validate_path<'a>(path: Option<&'a PathBuf>, name: &'a str, suffix: &str) -> &'a PathBuf {
    let path = path.expect(&format!("{} path is required", name));

    if !path.exists() {
        error!("Error: The specified path for {} does not exist", name);
        std::process::exit(1);
    }

    if !path.is_dir() {
        error!("Error: The specified path for {} is not a directory", name);
        exit(1);
    }

    let contains_files_with_suffix = fs::read_dir(path)
        .expect("Failed to read directory")
        // Filter out invalid entries
        .filter_map(|entry| entry.ok())
        .any(|entry| {
            if let Some(file_name) = entry.file_name().to_str() {
                file_name.contains(suffix)
            } else {
                false
            }
        });

    if !contains_files_with_suffix {
        error!(
            "Error: The directory for {} does not contain any files with the required extention/suffix '{}'",
            name, suffix
        );
        std::process::exit(1);
    }

    path
}

// Helper function
pub fn path_to_str(path: &PathBuf) -> &str {
    path.to_str().expect("Failed to convert PathBuf to &str")
}

// Helper function
pub fn check_paired_reads(directory: &PathBuf) -> bool {
    fs::read_dir(directory)
        .ok()
        .and_then(|entries| {
            entries
                .filter_map(|entry| entry.ok()?.file_name().to_str().map(String::from))
                .find(|name| name.contains("_1") || name.contains("_2"))
        })
        .is_some()
}

// Helper function
pub fn find_file_with_extension(directory: &PathBuf, base_name: &str) -> PathBuf {
    let fastq = directory.join(format!("{}.fastq", base_name));
    if fastq.exists() {
        fastq
    } else {
        directory.join(format!("{}.fastq.gz", base_name))
    }
}

// Helper function
pub fn get_binfiles(dir: &Path, extension: &str) -> io::Result<Vec<PathBuf>> {
    let mut files = Vec::new();

    for entry in fs::read_dir(dir)? {
        let entry = entry?;
        let path = entry.path();

        if path.is_file() {
            let filename = path
                .file_name()
                .unwrap_or_default()
                .to_str()
                .unwrap_or_default();
            if filename.contains("_all_seqs")
                || filename.contains("rep_seq")
                || filename.contains("combined")
            {
                continue;
            }

            if let Some(ext) = path.extension() {
                if ext == extension {
                    files.push(path);
                }
            }
        }
    }

    Ok(files)
}

// Helper function
pub fn get_sample_names(bindir: &Path, extension: &str) -> io::Result<HashMap<String, String>> {
    let mut bin_sample_map = HashMap::new();
    for entry in fs::read_dir(bindir)? {
        let path = entry?.path();
        if !path.is_file() || path.extension().and_then(|e| e.to_str()) != Some(extension) {
            continue;
        }
        // Open the file and read the first line
        let file_name = path
            .file_stem()
            .and_then(|name| name.to_str()) // Convert OsStr to &str
            .unwrap_or("unknown") // Fallback if conversion fails
            .to_string();
        let file = fs::File::open(&path)?;
        let mut reader = io::BufReader::new(file);
        let mut first_line = String::new();

        if reader.read_line(&mut first_line)? == 0 {
            continue;
        }

        let sample_id = first_line
            .trim()
            .trim_start_matches('>')
            .split('C')
            .next()
            .map(|s| s.to_string())
            .unwrap_or_default();

        bin_sample_map.insert(file_name, sample_id);
    }
    Ok(bin_sample_map)
}

// Helper function
pub fn splitbysampleid(
    bin: &PathBuf,
    bin_name: &str,
    binspecificdir: &PathBuf,
    format: &str,
) -> io::Result<()> {
    if !binspecificdir.exists() {
        fs::create_dir_all(&binspecificdir)?;
    }
    // Open the input file
    let reader = BufReader::new(File::open(&bin)?);

    // Create a HashMap to store writers for each sample ID
    let mut writers: HashMap<String, File> = HashMap::new();
    let mut current_sample_id = String::new();
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            current_sample_id = extract_sample_id(&line)?;
            ensure_writer(
                &current_sample_id,
                bin_name,
                binspecificdir,
                format,
                &mut writers,
            )?;
        }
        write_line_to_file(&current_sample_id, &line, &mut writers)?;
    }
    Ok(())
}

// Split bins by sample id and return the directory containing split bins.
pub fn split_bins_by_sample(
    parentdir: &Path,
    binfiles: &[PathBuf],
    format: &str,
    pool: &ThreadPool,
) -> io::Result<PathBuf> {
    let samplewisebinspath = parentdir.join("samplewisebins");
    if samplewisebinspath.exists() {
        fs::remove_dir_all(&samplewisebinspath)?;
    }
    fs::create_dir(&samplewisebinspath)?;

    pool.install(|| {
        binfiles
            .par_iter()
            .filter_map(|bin| bin.canonicalize().ok())
            .filter(|bin| {
                if bin.exists() {
                    true
                } else {
                    error!("Bin file does not exist: {:?}", bin);
                    false
                }
            })
            .for_each(|bin| {
                let bin_name = bin
                    .file_stem()
                    .and_then(|name| name.to_str())
                    .unwrap_or_default();
                if let Err(e) = splitbysampleid(&bin, bin_name, &samplewisebinspath, format) {
                    error!("Failed to split bin {:?}: {}", bin, e);
                }
            });
    });

    info!(
        "splitting bins by sample {:?} is completed",
        samplewisebinspath
    );
    Ok(samplewisebinspath)
}

// Helper function
pub fn extract_sample_id(line: &str) -> io::Result<String> {
    if let Some(idx) = line.find('C') {
        Ok(line[1..idx].to_string()) // Exclude the '>'
    } else {
        error!("Warning: Could not find 'C' in header: {}", line);
        Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "Invalid header format",
        ))
    }
}

// Helper function
pub fn ensure_writer(
    sample_id: &str,
    bin_name: &str,
    binspecificdir: &Path,
    format: &str,
    writers: &mut HashMap<String, File>,
) -> io::Result<()> {
    if !writers.contains_key(sample_id) {
        let output_filename = binspecificdir.join(format!("{}_{}.{}", bin_name, sample_id, format));
        let output_file = File::create(output_filename)?;
        writers.insert(sample_id.to_string(), output_file);
    }
    Ok(())
}

// Helper function
pub fn assign_members(
    component: &HashSet<String>,
    rep: &str,
    memberships_map: &Arc<Mutex<HashMap<String, String>>>,
) {
    if let Ok(mut map) = memberships_map.lock() {
        let rep_s = rep.to_string();
        for m in component.iter() {
            // member -> representative
            map.insert(m.clone(), rep_s.clone());
        }
    }
}

pub fn rep_members_from_member_rep(
    member_to_rep: &HashMap<String, String>,
) -> HashMap<String, Vec<String>> {
    let mut rep_to_members: HashMap<String, Vec<String>> = HashMap::new();

    for (member, rep) in member_to_rep {
        rep_to_members.entry(rep.clone()).or_default();
        if member != rep {
            rep_to_members
                .entry(rep.clone())
                .or_default()
                .push(member.clone());
        }
    }

    rep_to_members
}

pub fn write_membership_file(
    memberships_map: &HashMap<String, Vec<String>>,
    output_path: &Path,
) -> io::Result<()> {
    let output_file = File::create(output_path)?;
    let mut writer = io::BufWriter::new(output_file);
    let mut representatives: Vec<&String> = memberships_map.keys().collect();
    representatives.sort_unstable();

    writeln!(writer, "#representative\tmember_genomes")?;
    for rep in representatives {
        let mut members = memberships_map.get(rep).cloned().unwrap_or_default();
        members.sort_unstable();
        members.dedup();
        writeln!(writer, "{}\t{}", rep, members.join(","))?;
    }

    info!("Membership details are written to {:?}", output_path);
    Ok(())
}

// Helper function
pub fn write_line_to_file(
    sample_id: &str,
    line: &str,
    writers: &mut HashMap<String, File>,
) -> io::Result<()> {
    if let Some(writer) = writers.get_mut(sample_id) {
        writeln!(writer, "{}", line)?;
    }
    Ok(())
}

// Helper function
pub fn read_fasta(fasta_file: &str) -> io::Result<HashSet<String>> {
    let content = read_to_string(fasta_file)?;
    let mut scaffolds = HashSet::new();
    for line in content.lines() {
        if line.starts_with(">") {
            let scaffold_name = line
                .trim_start_matches(">")
                .split_whitespace()
                .next()
                .unwrap_or(line.trim_start_matches(">"));
            scaffolds.insert(scaffold_name.to_string());
        }
    }
    Ok(scaffolds)
}

// Helper function
pub fn get_output_binname(bin_fasta: &str) -> PathBuf {
    let path = Path::new(bin_fasta);

    let output_dir = path.parent().unwrap_or_else(|| Path::new("."));
    let filename = path
        .file_stem()
        .map(|stem| stem.to_str().unwrap_or("default"))
        .unwrap_or("default");
    output_dir.join(format!("{}.fasta", filename))
}
