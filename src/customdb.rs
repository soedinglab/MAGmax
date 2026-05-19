use std::collections::HashSet;
use std::fs;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::{Path, PathBuf};

use clap::Args;
use log::info;
use rayon::ThreadPoolBuilder;

use crate::{assess, utility};

#[derive(Args, Debug)]

pub struct CustomDbArgs {
    /// Directory containing fasta files of bins
    #[arg(
        short = 'b',
        long = "bindir",
        help = "Directory containing fasta files of bins"
    )]
    pub bindir: Option<PathBuf>,

    /// GTDB-Tk classification file
    #[arg(
        short = 'g',
        long = "gtdbtk",
        help = "GTDB-Tk classification summary file"
    )]
    pub gtdbtk: PathBuf,

    /// CheckM2 quality file
    #[arg(
        short = 'q',
        long = "qual",
        help = "Quality file produced by CheckM2 (quality_report.tsv)"
    )]
    pub qual: Option<PathBuf>,

    /// Select representatives based on high connectivity
    #[arg(
        long = "sensitive",
        help = "Select representatives based on high connectivity. Bin merging and reassembly steps are disabled"
    )]
    pub sensitive: bool,

    /// Average Nucleotide Identity for species-level clustering
    #[arg(
        long = "species-ani",
        default_value_t = 95.0,
        help = "ANI for clustering bins (%), as per GTDB-Tk criteria"
    )]
    pub species_ani: f32,

    /// Alignment fraction covered for species-level clustering
    #[arg(
        long = "species-alignedfrac",
        default_value_t = 50.0,
        help = "Minimum aligned fraction (%) for species-level clustering, as per GTDB-Tk criteria"
    )]
    pub species_alignedfrac: f32,

    /// Completeness cutoff
    #[arg(
        short = 'c',
        long = "completeness",
        default_value_t = 90.0,
        help = "Minimum completeness of bins (%)"
    )]
    pub completeness_cutoff: f32,

    /// Purity cutoff
    #[arg(
        short = 'p',
        long = "purity",
        default_value_t = 5.0,
        help = "Purity cutoff for custom database generation (%)"
    )]
    pub purity_cutoff: f32,

    /// Number of threads to use
    #[arg(
        short = 't',
        long = "threads",
        default_value_t = 8,
        help = "Number of threads to use"
    )]
    pub threads: usize,

    /// First split bins before processing
    #[arg(
        long = "split",
        help = "Split clusters into sample-wise bins before processing"
    )]
    pub split: bool,

    /// Bin file extension
    #[arg(
        short = 'f',
        long = "format",
        default_value = "fasta",
        help = "Bin file extension"
    )]
    format: String,

    /// ANI file
    #[arg(
        long = "anifile",
        help = "ANI file produced by skani using command: skani triangle <bindir> -E -o <anifile>"
    )]
    pub anifile: Option<PathBuf>,

    /// Directory of output
    #[arg(short = 'o', long = "outdir", help = "Directory of output")]
    pub output: Option<PathBuf>,
}

pub fn run(args: &CustomDbArgs) -> io::Result<()> {
    info!("Starting MAGmax custom database generation with parameters:");
    if let Some(bindir) = &args.bindir {
        info!("  Bins Directory: {:?}", bindir);
    }
    info!("  GTDB-Tk file: {:?}", args.gtdbtk);
    if args.sensitive {
        info!("  Mode: sensitive");
    } else {
        info!("  Mode: no-reassembly");
    }
    info!("  ANI Cutoff: {:.2}%", args.species_ani);
    info!(
        "  Aligned fraction cutoff: {:.1}%",
        args.species_alignedfrac
    );
    info!("  Completeness Cutoff: {:.1}%", args.completeness_cutoff);
    info!("  Purity Cutoff: {:.1}%", args.purity_cutoff);
    info!("  Threads: {}", args.threads);
    if let Some(qual) = &args.qual {
        info!("  CheckM2 quality file: {:?}", qual);
    }
    if let Some(anifile) = &args.anifile {
        info!("  ANI file: {:?}", anifile);
    }
    if let Some(output) = &args.output {
        info!("  Output Directory: {:?}", output);
    }

    let mut bindir = utility::validate_path(args.bindir.as_ref(), "bindir", &args.format).to_path_buf();
    
    if args.split {
        info!("  Split bins before processing: true");

        let parentdir = bindir
            .parent()
            .map(PathBuf::from)
            .unwrap_or_else(|| bindir.clone());
        let binfiles = utility::get_binfiles(&bindir, &args.format)?;
        if binfiles.is_empty() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "No bin files found. Please provide fasta files for customdb split.",
            ));
        }
        let pool = ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build()
            .expect("Failed to build thread pool");

        bindir = utility::split_bins_by_sample(&parentdir, &binfiles, &args.format, &pool)?;
    }

    let high_quality_bins = filter_bins_by_quality(args, Some(&bindir))?;
    info!(
        "Quality filtering retained {} bins",
        high_quality_bins.len()
    );

    let gtdb_bins = parse_gtdbtk_summary(
        &args.gtdbtk,
        args.species_ani,
        args.species_alignedfrac,
        &high_quality_bins,
    )?;
    let perfect_bins = gtdb_bins.perfect;
    let remaining_bins = gtdb_bins.remaining;

    info!(
        "GTDB-Tk summary parsed: {} perfect bins and {} remaining/unclassified bins",
        perfect_bins.len(),
        remaining_bins.len()
    );

    Ok(())
}

fn filter_bins_by_quality(
    args: &CustomDbArgs,
    bindir: Option<&PathBuf>,
) -> io::Result<HashSet<String>> {
    let checkm2_qualities = if let Some(qual_path) = &args.qual {
        if qual_path.is_file()
            && fs::metadata(qual_path)
                .map(|m| m.len() > 0)
                .unwrap_or(false)
        {
            qual_path.clone()
        } else {
            info!(
                "Provided quality file {:?} is missing or empty. Running CheckM2...",
                qual_path
            );
            run_checkm2_for_customdb(args, bindir)?
        }
    } else {
        run_checkm2_for_customdb(args, bindir)?
    };

    let bin_qualities = assess::parse_bins_quality(&checkm2_qualities).map_err(|e| {
        io::Error::new(
            e.kind(),
            format!(
                "Failed to parse CheckM2 quality file {:?}: {}",
                checkm2_qualities, e
            ),
        )
    })?;

    Ok(assess::filter_bins_quality(
        &bin_qualities,
        args.completeness_cutoff,
        args.purity_cutoff,
    ))
}

fn run_checkm2_for_customdb(args: &CustomDbArgs, bindir: Option<&PathBuf>) -> io::Result<PathBuf> {
    let bindir = bindir.ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            "--bindir is required when --qual is not provided or is missing/empty",
        )
    })?;
    let output_dir = args.output.clone().unwrap_or_else(|| {
        bindir
            .parent()
            .map(|parent| parent.join("customdb_results"))
            .unwrap_or_else(|| bindir.join("customdb_results"))
    });
    fs::create_dir_all(&output_dir)?;
    let checkm2_outputpath = output_dir.join("checkm2_results");

    assess::assess_bins(bindir, &checkm2_outputpath, args.threads, &args.format)
}

pub type PerfectGtdbBins = Vec<(String, String)>;

#[derive(Debug, Clone)]
pub struct GtdbBins {
    pub perfect: PerfectGtdbBins,
    pub remaining: Vec<String>,
}

fn normalize_af_cutoff(alignedfrac: f32) -> f32 {
    if alignedfrac > 1.0 {
        alignedfrac / 100.0
    } else {
        alignedfrac
    }
}

fn parse_gtdbtk_summary(
    summary_path: &Path,
    ani_species: f32,
    af_species: f32,
    high_quality_bins: &HashSet<String>,
) -> io::Result<GtdbBins> {
    const COL_GENOME: usize = 0;
    const COL_CLASSIF: usize = 1;
    const COL_CANI: usize = 5;
    const COL_CAF: usize = 6;

    let infile = File::open(summary_path)?;
    let reader = BufReader::new(infile);
    let mut perfect = Vec::new();
    let mut remaining = Vec::new();
    let af_species = normalize_af_cutoff(af_species);

    for (line_number, line) in reader.lines().enumerate() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if line_number == 0 && is_gtdbtk_header(&fields) {
            continue;
        }

        if fields.len() <= COL_CAF {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "GTDB-Tk summary row {} has {} columns; expected at least {}",
                    line_number + 1,
                    fields.len(),
                    COL_CAF + 1
                ),
            ));
        }

        let genome = fields[COL_GENOME].trim();
        let genome_name = strip_bin_extension(genome);
        if !high_quality_bins.contains(genome_name) {
            continue;
        }
        let classification = fields[COL_CLASSIF].trim();
        let closest_ani = parse_float(fields[COL_CANI]);
        let closest_af = normalize_af_cutoff(parse_float(fields[COL_CAF]));
        let species = assigned_species(classification);

        // checks if the bin meets GTDB-Tk species-level criteria:
        // ANI >= 95% and aligned fraction >= 50% with an assigned species name
        if closest_ani >= ani_species && closest_af >= af_species {
            if let Some(species) = species {
                perfect.push((genome_name.to_string(), species.to_string()));
                continue;
            }
        }
        remaining.push(genome_name.to_string());
    }

    Ok(GtdbBins { perfect, remaining })
}

fn strip_bin_extension(bin: &str) -> &str {
    for ext in [".fasta", ".faa", ".fna", ".fa", ".fas", ".ffn"] {
        if let Some(stripped) = bin.strip_suffix(ext) {
            return stripped;
        }
    }
    bin
}

fn is_gtdbtk_header(fields: &[&str]) -> bool {
    let first = fields.first().map(|field| field.trim()).unwrap_or_default();
    let second = fields.get(1).map(|field| field.trim()).unwrap_or_default();

    matches!(first, "user_genome" | "genome" | "genome_id") || second == "classification"
}

fn parse_float(value: &str) -> f32 {
    value.trim().parse::<f32>().unwrap_or(-1.0)
}

fn assigned_species(classification: &str) -> Option<&str> {
    let species = classification.rsplit("s__").next()?.trim();
    if species.is_empty() || species == classification {
        None
    } else {
        Some(species)
    }
}
