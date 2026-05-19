use csv::ReaderBuilder;
use log::error;
use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{self};
use std::path::{Path, PathBuf};
use std::process::{Command as ProcessCommand, Stdio};

#[derive(Clone, Debug)]
pub struct BinQuality {
    pub completeness: f32,
    pub contamination: f32,
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
            Ok(_) => {}
            Err(e) => {
                error!(
                    "Error: Failed to execute CheckM2 command - {}. 
                    Check if CheckM2 is executable currently",
                    e
                );
            }
        }
    }
    Ok(checkm2_qualities)
}

/// Parse CheckM2 result
pub fn parse_bins_quality(checkm2_qualities: &PathBuf) -> io::Result<HashMap<String, BinQuality>> {
    let _ = File::open(checkm2_qualities).map_err(|e| {
        io::Error::new(
            e.kind(),
            format!("Failed to open checkm2 quality file for bin: {:?}", e),
        )
    })?;

    let mut rdr = ReaderBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .from_reader(File::open(checkm2_qualities)?);
    let mut bin_qualities: HashMap<String, BinQuality> = HashMap::new();

    for result in rdr.records() {
        let record = result?; // Get the record from the CSV
        if record.len() < 3 {
            error!("Skipping invalid record: {:?}", record);
            continue; // Skip records that do not have enough columns
        }
        let mut bin_id: String = record[0].to_string();
        for ext in [".fasta", ".faa", ".fna", ".fa", ".fas", ".ffn"] {
            if let Some(stripped) = bin_id.strip_suffix(ext) {
                bin_id = stripped.to_string();
                break;
            }
        }
        let completeness: f32 = record[1].parse().unwrap_or(0.0);
        let contamination: f32 = record[2].parse().unwrap_or(0.0);
        bin_qualities.insert(
            bin_id,
            BinQuality {
                completeness,
                contamination,
            },
        );
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
) -> Option<String> {
    let comp_binqualities: HashMap<String, BinQuality> = bin_qualities
        .iter()
        .filter(|(bin, _)| comp.contains(bin.as_str()))
        .map(|(bin, q)| (bin.clone(), q.clone()))
        .collect();

    if let Some((binname, _)) = comp_binqualities
        .iter()
        .filter(|(_, q)| q.completeness > 90.0)
        .max_by(|(bin1, q1), (bin2, q2)| {
            // select the best bin by quality score
            let score1 = q1.completeness - (5.0 * q1.contamination);
            let score2 = q2.completeness - (5.0 * q2.contamination);

            score1
                .total_cmp(&score2)
                .then_with(|| q2.contamination.total_cmp(&q1.contamination))
                .then_with(|| bin2.cmp(bin1).reverse()) // Tie-breaker based on contamination
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
