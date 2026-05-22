use csv::ReaderBuilder;
use std::borrow::Cow;
use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};
use std::fs::{self, File};
use std::io::{self};
use std::process::{Command as ProcessCommand, Stdio};
use log::error;

#[derive(Clone)]
#[derive(Debug)]
pub struct BinQuality {
    pub completeness: f32,
    pub contamination: f32,
}

impl BinQuality {
    pub fn score(&self) -> f32 {
        quality_score(self.completeness, self.contamination)
    }
}

pub fn quality_score(completeness: f32, contamination: f32) -> f32 {
    completeness - 5.0 * contamination
}

/// Run CheckM2 to obtain completeness and contamination of input bins
pub fn assess_bins(
    bindir: &PathBuf,
    bincheckm2dir: &PathBuf,
    threads: usize,
    format: &str,
) -> Result<PathBuf, io::Error> {

    let checkm2_qualities = Path::new(bincheckm2dir).join("quality_report.tsv");

    if !checkm2_qualities.exists() {
        let mut output = ProcessCommand::new("checkm2");
        output
        .arg("predict")
        .arg("-i")
        .arg(bindir)
        .arg("-o")
        .arg(bincheckm2dir)
        .arg("-t")
        .arg(threads.to_string())
        .arg("-x")
        .arg(format)
        .arg("--force")
        .stdout(Stdio::null())
        .stderr(Stdio::null());
    
        match output.status() {
            Ok(status) if status.success() => {}
            Ok(status) => {
                return Err(io::Error::new(
                    io::ErrorKind::Other,
                    format!("CheckM2 exited with non-zero status: {}", status),
                ));
            }
            Err(e) => {
                return Err(io::Error::new(
                    e.kind(),
                    format!("Failed to execute CheckM2: {}. Check if CheckM2 is in PATH", e),
                ));
            }
        }
    }
    Ok(checkm2_qualities)
}

/// Parse CheckM2 result
pub fn parse_bins_quality(
    checkm2_qualities: &PathBuf,
) -> io::Result<HashMap<String, BinQuality>> {

    let file = File::open(checkm2_qualities).map_err(|e| {
        io::Error::new(
            e.kind(),
            format!("Failed to open checkm2 quality file for bin: {:?}", e),
        )
    })?;

    let mut rdr = ReaderBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .from_reader(file);
    let mut bin_qualities: HashMap<String, BinQuality> = HashMap::new();

    for result in rdr.records() {
        let record = result?; // Get the record from the CSV
        if record.len() < 3 {
            error!("Skipping invalid record: {:?}", record);
            continue; // Skip records that do not have enough columns
        }
        let raw: &str = &record[0];
        let bin_id: Cow<str> = [".fasta", ".faa", ".fna", ".fa", ".fas", ".ffn"]
            .iter()
            .find_map(|ext| raw.strip_suffix(ext))
            .map(Cow::Borrowed)
            .unwrap_or(Cow::Borrowed(raw));
        let completeness: f32 = record[1].parse().unwrap_or(0.0);
        let contamination: f32 = record[2].parse().unwrap_or(0.0);
        bin_qualities.insert(bin_id.into_owned(),
        BinQuality{ completeness, contamination });
    }
    Ok(bin_qualities)
}

/// Filter bins by CheckM2 completeness and contamination thresholds.
pub fn filter_bins_quality(
    bin_qualities: &HashMap<String, BinQuality>,
    completeness_cutoff: f32,
    contamination_cutoff: f32,
) -> HashMap<String, BinQuality> {
    bin_qualities
        .iter()
        .filter(|(_, q)| {
            q.completeness >= completeness_cutoff && q.contamination <= contamination_cutoff
        })
        .map(|(bin, q)| (bin.clone(), q.clone()))
        .collect()
}


/// Select a high-quality bin if it exist in the cluster
pub fn check_high_quality_bin(
    comp: &HashSet<String>,
    bin_qualities: &HashMap<String, BinQuality>,
    bindir: &PathBuf,
    resultdir: &PathBuf,
    format: &str,
    completeness_cutoff: f32,
) -> Option<String> {
    let hq_threshold = completeness_cutoff.max(90.0);

    if let Some((binname, _)) = comp
        .iter()
        .filter_map(|bin| bin_qualities.get(bin).map(|q| (bin, q)))
        .filter(|(_, q)| q.completeness > hq_threshold)
        .max_by(|(_, q1), (_, q2)| {
            q1.score()
                .total_cmp(&q2.score())
                .then_with(|| q2.contamination.total_cmp(&q1.contamination))
        })

    {
        let bin_path = bindir.join(format!("{}.{}", binname, format));
        let final_path = resultdir.join(format!("{}.{}", binname, format));

        if let Err(_) = fs::copy(&bin_path, &final_path) {
            return None;
        }
        return Some(binname.clone());
    }
    None
}
