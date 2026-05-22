use std::collections::{HashMap, HashSet};
use std::fs;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use clap::Args;
use log::{debug, info, warn};
use rayon::ThreadPoolBuilder;

use crate::{assess, merge, mwids, reassemble, utility};

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

    /// File listing isolate genomes present in the input bin directory
    #[arg(
        long = "isolate-genomes",
        help = "File listing isolate genomes in the input bins; these are prioritized as species representatives"
    )]
    pub isolate_genomes: Option<PathBuf>,

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
    if let Some(isolate_genomes) = &args.isolate_genomes {
        info!("  Isolate genome list: {:?}", isolate_genomes);
    }
    if let Some(anifile) = &args.anifile {
        info!("  ANI file: {:?}", anifile);
    }
    if let Some(output) = &args.output {
        info!("  Output Directory: {:?}", output);
    }

    let mut bindir =
        utility::validate_path(args.bindir.as_ref(), "bindir", &args.format).to_path_buf();

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

    let bin_qualities = filter_bins_by_quality(args, Some(&bindir))?;
    let quality_filtered_bins: HashSet<String> = bin_qualities.keys().cloned().collect();
    info!(
        "Quality filtering retained {} bins",
        quality_filtered_bins.len()
    );

    let isolate_genomes = read_isolate_genomes(args.isolate_genomes.as_deref())?;
    let ignored_isolates = isolate_genomes
        .len()
        .saturating_sub(isolate_genomes.intersection(&quality_filtered_bins).count());

    if ignored_isolates > 0 {
        warn!(
            "{} isolate genomes from the input are not quality passed and will be ignored",
            ignored_isolates
        );
    }

    let gtdb_bins = parse_gtdbtk_summary(
        &args.gtdbtk,
        args.species_ani,
        args.species_alignedfrac,
        &quality_filtered_bins,
    )?;

    let perfect_bins = gtdb_bins.perfect;
    let remaining_bins = gtdb_bins.remaining;
    let sp_aniradius = gtdb_bins.sp_aniradius;

    info!(
        "GTDB-Tk summary parsed: {} perfect bins and {} remaining/unclassified bins",
        perfect_bins.len(),
        remaining_bins.len()
    );

    let output_dir = customdb_output_dir(args, &bindir);
    fs::create_dir_all(&output_dir)?;

    let mut memberships_map =
        group_perfect_bins_by_species(&perfect_bins, &bin_qualities, &isolate_genomes);
    write_gtdbtk_species_representatives(
        &memberships_map,
        &sp_aniradius,
        &output_dir.join("gtdbtk_species_representatives.tsv"),
    )?;
    let mut representative_bins: HashSet<String> = memberships_map.keys().cloned().collect();

    let remaining_bin_set: HashSet<String> = remaining_bins.into_iter().collect();
    let remaining_qualities: HashMap<String, assess::BinQuality> = bin_qualities
        .iter()
        .filter(|(bin, _)| remaining_bin_set.contains(*bin))
        .map(|(bin, quality)| (bin.clone(), quality.clone()))
        .collect();

    if !remaining_qualities.is_empty() {
        let (remaining_reps, remaining_memberships) = process_remaining_bins(
            args,
            &bindir,
            &output_dir,
            &representative_bins,
            &remaining_qualities,
            &isolate_genomes,
        )?;
        write_remaining_species_ani_report(
            args,
            &bindir,
            &remaining_reps,
            &mut memberships_map,
            &representative_bins,
            &sp_aniradius,
            &output_dir,
        )?;
        representative_bins.extend(remaining_reps);
        merge_memberships(&mut memberships_map, remaining_memberships);
    }

    // copy_representative_bins(&representative_bins, &bindir, &output_dir, &args.format)?;
    utility::write_membership_file(&memberships_map, &output_dir.join("memberships.tsv"))?;
    write_quality_file(&representative_bins, &bin_qualities, &output_dir)?;

    Ok(())
}

fn write_remaining_species_ani_report(
    args: &CustomDbArgs,
    bindir: &PathBuf,
    remaining_reps: &HashSet<String>,
    memberships_map: &mut HashMap<String, Vec<String>>,
    representative_bins: &HashSet<String>,
    sp_aniradius: &HashMap<String, (f32, String)>,
    output_dir: &Path,
) -> io::Result<()> {
    if remaining_reps.is_empty() || representative_bins.is_empty() {
        return Ok(());
    }

    // `ani_details` is keyed by internal node IDs from the remaining-bin graph.
    // For remaining-vs-existing-cluster checks, use the skani TSV keyed by bin names.
    let (ani_by_name, afr_by_name, afq_by_name) = load_ani_by_name(args, bindir)?;

    let report_path =
        output_dir.join("unclassified_clusterrepresentatives_gtdbtkspecies_ani_connections.tsv");
    let report_file = File::create(&report_path)?;
    let mut writer = BufWriter::new(report_file);
    let mut reps: Vec<&String> = remaining_reps.iter().collect();
    reps.sort_unstable();
    let mut cluster_reps: Vec<&String> = representative_bins.iter().collect();
    cluster_reps.sort_unstable();

    writeln!(
        writer,
        "#unclassified_cluster_representative\tgtdbtk_species_representative\tANI\tspecies_ANI_radius\tall_members_share_ANI_above_radius"
    )?;

    for remaining_rep in reps {
        for cluster_rep in &cluster_reps {
            if remaining_rep == *cluster_rep {
                warn!(
                    "Unclassified representative {} is identical to GTDB-Tk species cluster representative {}. Check carefully!",
                    remaining_rep, cluster_rep
                );
                continue;
            }

            let mut cluster_members = Vec::new();
            cluster_members.push((*cluster_rep).clone());
            if let Some(members) = memberships_map.get(*cluster_rep) {
                cluster_members.extend(members.iter().cloned());
            }
            cluster_members.sort_unstable();

            let mut all_members_above_species_ani = true;

            let Some(rep_ani) = lookup_ani_by_name(remaining_rep, cluster_rep, &ani_by_name) else {
                continue;
            };
            debug!(
                "species ANI radius: {}",
                sp_aniradius.get(*cluster_rep).map(|(ani_radius, _)| ani_radius.to_string())
                    .unwrap_or_else(|| "unknown".to_string())
            );
            if rep_ani < sp_aniradius.get(*cluster_rep).map(|(ani_radius, _)| *ani_radius).unwrap_or(0.0) {
                continue;
            }
            let key = ordered_name_pair(remaining_rep.clone(), (*cluster_rep).clone());
            if afr_by_name.get(&key).copied().unwrap_or(0.0) < args.species_alignedfrac {
                continue;
            }
            if afq_by_name.get(&key).copied().unwrap_or(0.0) < args.species_alignedfrac {
                continue;
            }
            for member in &cluster_members {
                if remaining_rep == member {
                    warn!(
                    "Unclassified representative {} is identical to GTDB-Tk species cluster member {}. Check carefully!",
                    remaining_rep, member
                    );
                    continue;
                }

                let Some(ani) = lookup_ani_by_name(remaining_rep, member, &ani_by_name) else {
                    all_members_above_species_ani = false;
                    break;
                };

                if ani < sp_aniradius.get(member).map(|(ani_radius, _)| *ani_radius).unwrap_or(0.0) {
                    all_members_above_species_ani = false;
                    break;
                }
                let key = ordered_name_pair(remaining_rep.clone(), member.clone());
                if afr_by_name.get(&key).copied().unwrap_or(0.0) < args.species_alignedfrac {
                    all_members_above_species_ani = false;
                    break;
                }
                if afq_by_name.get(&key).copied().unwrap_or(0.0) < args.species_alignedfrac {
                    all_members_above_species_ani = false;
                    break;
                }
            }

            writeln!(
                writer,
                "{}\t{}\t{:4}\t{:4}\t{}",
                remaining_rep,
                cluster_rep,
                rep_ani,
                sp_aniradius.get(*cluster_rep).map(|(ani_radius, _)| *ani_radius).unwrap_or(0.0),
                all_members_above_species_ani,
            )?;
        }
    }

    info!(
        "Remaining representative species-ANI report is written to {:?}",
        report_path
    );
    Ok(())
}

fn load_ani_by_name(
    args: &CustomDbArgs,
    bindir: &PathBuf,
) -> io::Result<(
    HashMap<(String, String), f32>,
    HashMap<(String, String), f32>,
    HashMap<(String, String), f32>,
)> {
    let ani_file = args
        .anifile
        .clone()
        .filter(|path| path.exists())
        .unwrap_or_else(|| bindir.join("ani_edges"));

    let infile = File::open(&ani_file)?;
    let reader = BufReader::new(infile);
    let mut ani_by_name = HashMap::new();
    let mut afr_by_name = HashMap::new();
    let mut afq_by_name = HashMap::new();
    for line in reader.lines().skip(1) {
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 5 {
            continue;
        }

        let bin1 = ani_bin_name(fields[0]);
        let bin2 = ani_bin_name(fields[1]);
        let ani = fields[2].parse::<f32>().unwrap_or(0.0);
        let afr = fields[3].parse::<f32>().unwrap_or(0.0);
        let afq = fields[4].parse::<f32>().unwrap_or(0.0);
        let key = ordered_name_pair(bin1, bin2);
        ani_by_name.insert(key.clone(), ani);
        afr_by_name.insert(key.clone(), afr);
        afq_by_name.insert(key, afq);
    }

    Ok((ani_by_name, afr_by_name, afq_by_name))
}

fn ani_bin_name(value: &str) -> String {
    Path::new(value)
        .file_stem()
        .map(|name| name.to_string_lossy().into_owned())
        .unwrap_or_else(|| strip_bin_extension(value).to_string())
}

fn lookup_ani_by_name(
    bin1: &str,
    bin2: &str,
    ani_by_name: &HashMap<(String, String), f32>,
) -> Option<f32> {
    ani_by_name
        .get(&ordered_name_pair(bin1.to_string(), bin2.to_string()))
        .copied()
}

fn get_remaining_ani_file(
    args: &CustomDbArgs,
    bindir: &PathBuf,
    output_dir: &PathBuf,
    representative_bins: &HashSet<String>,
    remaining_qualities: &HashMap<String, assess::BinQuality>,
) -> io::Result<PathBuf> {
    if let Some(anifile) = &args.anifile {
        return Ok(anifile.clone());
    }

    let ani_output = output_dir.join("ani_edges");
    if ani_output.exists() {
        fs::remove_file(&ani_output).ok();
    }

    let mut input_bin_paths = Vec::new();
    for bin in representative_bins {
        let bin_path = bindir.join(format!("{}.{}", bin, args.format));
        if !bin_path.exists() {
            return Err(io::Error::new(
                io::ErrorKind::NotFound,
                format!("Perfect bin file not found: {:?}", bin_path),
            ));
        }
        input_bin_paths.push(bin_path.to_string_lossy().into_owned());
    }

    for bin in remaining_qualities.keys() {
        let bin_path = bindir.join(format!("{}.{}", bin, args.format));
        if !bin_path.exists() {
            return Err(io::Error::new(
                io::ErrorKind::NotFound,
                format!("Remaining bin file not found: {:?}", bin_path),
            ));
        }
        input_bin_paths.push(bin_path.to_string_lossy().into_owned());
    }

    if input_bin_paths.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "No input bin files available to generate ANI file",
        ));
    }

    info!(
        "Calculating ANI among perfect bins and remaining bins using skani: {:?}",
        ani_output
    );
    merge::get_ani(input_bin_paths, &ani_output, args.threads)?;
    Ok(ani_output)
}

fn ordered_name_pair(bin1: String, bin2: String) -> (String, String) {
    if bin1 <= bin2 {
        (bin1, bin2)
    } else {
        (bin2, bin1)
    }
}

fn process_remaining_bins(
    args: &CustomDbArgs,
    bindir: &PathBuf,
    output_dir: &PathBuf,
    representative_bins: &HashSet<String>,
    remaining_qualities: &HashMap<String, assess::BinQuality>,
    isolate_genomes: &HashSet<String>,
) -> io::Result<(HashSet<String>, HashMap<String, Vec<String>>)> {
    let ani_file = get_remaining_ani_file(args, bindir, output_dir, representative_bins, remaining_qualities)?;
    let (graph, ani_data, id_to_name, id_to_node) = merge::calc_ani(
        bindir,
        remaining_qualities,
        &args.format,
        Some(ani_file),
        args.species_ani,
        args.completeness_cutoff,
        args.purity_cutoff,
        args.species_alignedfrac,
        true,
        args.threads,
    )?;

    if args.sensitive {
        info!("Selecting remaining-bin representatives based on high connectivity");
        let ani_map: HashMap<(u32, u32), f32> = ani_data.iter()
            .map(|(&k, v)| (k, v.ani))
            .collect();

        return Ok(mwids::select_highconnectivity_bins_with_memberships(
            &graph,
            &ani_map,
            args.species_ani,
            &id_to_name,
            &id_to_node,
            remaining_qualities,
            isolate_genomes,
        ));
    }

    info!("Selecting remaining-bin representatives without reassembly");
    let connected_bins = merge::get_connected_samples(
        &graph,
        &ani_data,
        args.species_ani,
        &id_to_name,
        args.species_alignedfrac,
    );

    let mut representative_bins = HashSet::new();
    let mut member_to_rep = HashMap::new();

    for component in connected_bins {
        if let Some(rep) = select_no_reassembly_rep(
            &component,
            remaining_qualities,
            args.completeness_cutoff,
            isolate_genomes,
        ) {
            representative_bins.insert(rep.clone());
            for member in component {
                member_to_rep.insert(member, rep.clone());
            }
        }
    }

    Ok((
        representative_bins,
        utility::rep_members_from_member_rep(&member_to_rep),
    ))
}

fn select_no_reassembly_rep(
    component: &HashSet<String>,
    bin_qualities: &HashMap<String, assess::BinQuality>,
    completeness_cutoff: f32,
    isolate_genomes: &HashSet<String>,
) -> Option<String> {
    if component.len() == 1 {
        let binname = component.iter().next()?;
        let quality = bin_qualities.get(binname)?;
        if quality.completeness >= completeness_cutoff {
            return Some(binname.clone());
        }
        return None;
    }

    let isolate_bins = bins_in_isolate_list(component.iter(), isolate_genomes);
    if let Some(rep) = select_best_quality_bin(&isolate_bins, bin_qualities) {
        return Some(rep);
    }

    let high_quality_bins = component
        .iter()
        .filter(|bin| {
            bin_qualities
                .get(*bin)
                .map(|quality| quality.completeness > 90.0)
                .unwrap_or(false)
        })
        .cloned()
        .collect::<Vec<_>>();

    if let Some(rep) = select_best_quality_bin(&high_quality_bins, bin_qualities) {
        return Some(rep);
    }

    reassemble::find_bestqualitybin(component, bin_qualities, completeness_cutoff)
        .map(|(bin_name, _, _)| bin_name)
}

fn merge_memberships(
    memberships_map: &mut HashMap<String, Vec<String>>,
    additional_memberships: HashMap<String, Vec<String>>,
) {
    for (rep, members) in additional_memberships {
        memberships_map.entry(rep).or_default().extend(members);
    }
}

// fn copy_representative_bins(
//     representative_bins: &HashSet<String>,
//     bindir: &Path,
//     output_dir: &Path,
//     format: &str,
// ) -> io::Result<()> {
//     for bin in representative_bins {
//         let bin_path = bindir.join(format!("{}.{}", bin, format));
//         let final_path = output_dir.join(format!("{}.{}", bin, format));
//         if !final_path.exists() {
//             fs::copy(&bin_path, &final_path)?;
//         }
//     }
//     Ok(())
// }

fn write_gtdbtk_species_representatives(
    memberships_map: &HashMap<String, Vec<String>>,
    sp_aniradius: &HashMap<String, (f32, String)>,
    output_path: &Path,
) -> io::Result<()> {
    let output_file = File::create(output_path)?;
    let mut writer = BufWriter::new(output_file);
    let mut representatives: Vec<&String> = memberships_map.keys().collect();
    representatives.sort_unstable();

    writeln!(writer, "#gtdbtk_species_representative\tspecies_name")?;
    for rep in representatives {
        let sp_name = sp_aniradius
            .get(rep)
            .map(|(_, classification)| classification.clone())
            .unwrap_or_else(|| "unknown".to_string());
        writeln!(writer, "{}\t{}", rep, sp_name)?;
    }

    info!(
        "GTDB-Tk species representatives are written to {:?}",
        output_path
    );
    Ok(())
}

fn write_quality_file(
    representative_bins: &HashSet<String>,
    bin_qualities: &HashMap<String, assess::BinQuality>,
    output_dir: &Path,
) -> io::Result<()> {
    let output_file_path = output_dir.join("bins_checkm2_qualities.tsv");
    let output_file = File::create(&output_file_path)?;
    let mut writer = BufWriter::new(output_file);
    let mut representatives: Vec<&String> = representative_bins.iter().collect();
    representatives.sort_unstable();

    writeln!(writer, "#Bin\tCompleteness\tContamination")?;
    for bin in representatives {
        if let Some(quality) = bin_qualities.get(bin) {
            writeln!(
                writer,
                "{}\t{}\t{}",
                bin, quality.completeness, quality.contamination
            )?;
        }
    }
    info!(
        "Quality values of bins are written to {:?}",
        output_file_path
    );
    Ok(())
}

fn customdb_output_dir(args: &CustomDbArgs, bindir: &PathBuf) -> PathBuf {
    args.output.clone().unwrap_or_else(|| {
        bindir
            .parent()
            .map(|parent| parent.join("specieslevel_customdb"))
            .unwrap_or_else(|| bindir.join("specieslevel_customdb"))
    })
}

fn group_bins_by_species(perfect_bins: &PerfectGtdbBins) -> HashMap<String, Vec<String>> {
    let mut species_bins: HashMap<String, Vec<String>> = HashMap::new();

    for (bin, species) in perfect_bins {
        species_bins
            .entry(species.clone())
            .or_default()
            .push(bin.clone());
    }

    for bins in species_bins.values_mut() {
        bins.sort_unstable();
    }

    species_bins
}

fn group_perfect_bins_by_species(
    perfect_bins: &PerfectGtdbBins,
    bin_qualities: &HashMap<String, assess::BinQuality>,
    isolate_genomes: &HashSet<String>,
) -> HashMap<String, Vec<String>> {
    let mut species_bins = group_bins_by_species(perfect_bins);

    let mut memberships_map = HashMap::new();
    for bins in species_bins.values_mut() {
        bins.sort_unstable();
        let Some(rep) = select_isolate_priority_rep(bins, bin_qualities, isolate_genomes) else {
            continue;
        };

        let members = bins
            .iter()
            .filter(|bin| *bin != &rep)
            .cloned()
            .collect::<Vec<_>>();
        memberships_map.insert(rep, members);
    }

    memberships_map
}

fn select_isolate_priority_rep(
    bins: &[String],
    bin_qualities: &HashMap<String, assess::BinQuality>,
    isolate_genomes: &HashSet<String>,
) -> Option<String> {
    let isolate_bins = bins_in_isolate_list(bins.iter(), isolate_genomes);
    let best_isolate = select_best_quality_bin(&isolate_bins, bin_qualities);
    let best_overall = select_best_quality_bin(bins, bin_qualities);

    // select best purity bin among isolates if exist
    // otherwise select best purity bin among all bins in the species cluster
    match (&best_isolate, &best_overall) {
        (Some(isolate_bin), Some(overall_bin)) => {
            let isolate_cont = bin_qualities[isolate_bin].contamination;
            let overall_cont = bin_qualities[overall_bin].contamination;

            if overall_cont < isolate_cont {
                best_overall
            } else {
                best_isolate
            }
        }
        (Some(_), None) => best_isolate,
        (None, Some(_)) => best_overall,
        (None, None) => None,
    }
}

fn bins_in_isolate_list<'a, I>(bins: I, isolate_genomes: &HashSet<String>) -> Vec<String>
where
    I: IntoIterator<Item = &'a String>,
{
    bins.into_iter()
        .filter(|bin| isolate_genomes.contains(*bin))
        .cloned()
        .collect()
}

fn select_best_quality_bin(
    bins: &[String],
    bin_qualities: &HashMap<String, assess::BinQuality>,
) -> Option<String> {
    bins.iter()
        .filter_map(|bin| bin_qualities.get(bin).map(|quality| (bin, quality)))
        .max_by(|(bin1, q1), (bin2, q2)| {
            let score1 = q1.completeness - (5.0 * q1.contamination);
            let score2 = q2.completeness - (5.0 * q2.contamination);

            score1
                .total_cmp(&score2)
                .then_with(|| q2.contamination.total_cmp(&q1.contamination))
                .then_with(|| bin2.cmp(bin1).reverse())
        })
        .map(|(bin, _)| bin.clone())
}

fn read_isolate_genomes(isolate_path: Option<&Path>) -> io::Result<HashSet<String>> {
    let Some(isolate_path) = isolate_path else {
        return Ok(HashSet::new());
    };

    let infile = File::open(isolate_path)?;
    let reader = BufReader::new(infile);
    let mut isolate_genomes = HashSet::new();

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        if let Some(genome) = trimmed.split_whitespace().next() {
            isolate_genomes.insert(normalize_bin_name(genome));
        }
    }

    Ok(isolate_genomes)
}

fn normalize_bin_name(value: &str) -> String {
    let value = value.trim_end_matches(['/', '\\']);
    let basename = value.rsplit(['/', '\\']).next().unwrap_or(value);
    strip_bin_extension(basename).to_string()
}

fn filter_bins_by_quality(
    args: &CustomDbArgs,
    bindir: Option<&PathBuf>,
) -> io::Result<HashMap<String, assess::BinQuality>> {
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
    let output_dir = customdb_output_dir(args, bindir);
    fs::create_dir_all(&output_dir)?;
    let checkm2_outputpath = output_dir.join("checkm2_results");

    assess::assess_bins(bindir, &checkm2_outputpath, args.threads, &args.format)
}

pub type PerfectGtdbBins = HashMap<String, String>;

#[derive(Debug, Clone)]
pub struct GtdbBins {
    pub perfect: PerfectGtdbBins,
    pub remaining: Vec<String>,
    pub sp_aniradius: HashMap<String, (f32, String)>,
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
    quality_filtered_bins: &HashSet<String>,
) -> io::Result<GtdbBins> {
    const COL_GENOME: usize = 0;
    const COL_CLASSIF: usize = 1;
    const COL_ANIRAD: usize = 3;
    const COL_CANI: usize = 5;
    const COL_CAF: usize = 6;

    let infile = File::open(summary_path)?;
    let reader = BufReader::new(infile);
    let mut perfect = HashMap::new();
    let mut remaining = Vec::new();
    let mut sp_aniradius = HashMap::new();
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
        if !quality_filtered_bins.contains(genome_name) {
            continue;
        }
        let classification = fields[COL_CLASSIF].trim();
        let closest_ani = parse_float(fields[COL_CANI]);
        let closest_af = normalize_af_cutoff(parse_float(fields[COL_CAF]));
        let species = assigned_species(classification);
        let species_ani_radius = parse_float(fields[COL_ANIRAD]);
        let species_ani_radius = if species_ani_radius >= 0.0 {
            species_ani_radius
        } else {
            ani_species
        };

        if let Some(_) = species {
            sp_aniradius.insert(genome_name.to_string(), (species_ani_radius, classification.to_string()));
        }
        // checks if the bin meets GTDB-Tk species-level criteria:
        // ANI >= species cutoff, ANI >= GTDB-Tk species ANI radius, and aligned fraction >= cutoff
        if closest_ani >= species_ani_radius
            && closest_af >= af_species
        {
            if let Some(species) = species {
                perfect.insert(genome_name.to_string(), species.to_string());
                continue;
            }
        }
        remaining.push(genome_name.to_string());
    }

    Ok(GtdbBins {
        perfect,
        remaining,
        sp_aniradius,
    })
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
