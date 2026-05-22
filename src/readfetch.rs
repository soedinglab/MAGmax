use std::collections::HashSet;
use std::process::{Command, Stdio};
use std::path::PathBuf;
use std::fs::{File, OpenOptions};
use std::io::{self, BufRead, Write};
use std::thread;
use flate2::read::GzDecoder;
use crate::utility;
use log::error;

/// Read fastq file and collect reads mapped to contigs of merged bin
pub fn fetch_fastqreads(
    enriched_scaffolds: &HashSet<String>,
    mapids: &str,
    fastq_files: Vec<String>,
    outputbin: PathBuf,
    is_paired: bool
) -> Result<(), Box<dyn std::error::Error>> {

    let base = utility::get_output_binname(
            outputbin.to_str().expect("Invalid UTF-8 in outputbin path"))
        .to_str()
        .expect("Invalid UTF-8 in output bin name")
        .replace(".fasta", "");

    let output_fastq: Vec<String> = if is_paired {
        vec![format!("{}_1.fastq", base), format!("{}_2.fastq", base)]
    } else {
        vec![format!("{}.fastq", base)]
    };

    if let Err(e) = write_selected_reads(
        fastq_files,
        enriched_scaffolds,
        mapids,
        &output_fastq,
        is_paired
    ) {
        error!(" {}", e);
    };

    let file_metadata = std::fs::metadata(&output_fastq[0]).map_err(|e| {
        error!("Error accessing {:?}: {}", output_fastq, e);
        Box::new(std::io::Error::new(
            std::io::ErrorKind::Other,
            format!("Failed to access {:?}", output_fastq),
        )) as Box<dyn std::error::Error>
    })?;

    if file_metadata.len() == 0 {
        error!("Error: The output file {:?} is empty", output_fastq[0]);
        return Err(Box::new(std::io::Error::new(
            std::io::ErrorKind::Other,
            format!("No reads were written to {:?}", output_fastq[0]),
        )));
    }
    Ok(())
}

// Write mapped reads to fastq files
fn write_selected_reads(
    fastq_files: Vec<String>,
    enriched_scaffolds: &HashSet<String>,
    mapid_file: &str,
    output_fastq: &[String],
    is_paired: bool,
) -> Result<(), io::Error> {
    let mfile = File::open(mapid_file)?;
    let mapid_reader = io::BufReader::new(mfile);
    let readid_file = format!("{}_readids", output_fastq[0].replace(".fastq", ""));
    let mut idfile = io::BufWriter::new(
        File::create(&readid_file)?);
    for line in mapid_reader.lines() {
        let line = line?;
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 2 {
            error!("map id file is not in two column format");
            continue;
        }

        // As of now, it only works for readid format: @SRR3961047.1
        let read_id = parts[0];
        let scaffold_id = parts[1];
        // Check if scaffold_id exists in the enriched/combined set
        if enriched_scaffolds.contains(scaffold_id) {
            writeln!(idfile, "{}", read_id)?;
        }
    }
    idfile.flush()?;

    let mut outfile1 = OpenOptions::new()
        .create(true)
        .append(true)
        .open(&output_fastq[0])?;

    let mut outfile2 = if is_paired {
        Some(OpenOptions::new()
            .create(true)
            .append(true)
            .open(&output_fastq[1])?)
    } else {
        None
    };

    if which::which("seqtk").is_err() {
        return Err(io::Error::new(io::ErrorKind::NotFound, "`seqtk` not found in PATH"));
    }
    let process_seqtk =
        |fastq: &String, outfile: &mut File| -> Result<(), io::Error> {
        if fastq.ends_with(".gz") {
            // Decompress in a background thread piped to seqtk stdin to avoid
            // deadlock between the stdin write and stdout read.
            let mut child = Command::new("seqtk")
                .arg("subseq")
                .arg("-")
                .arg(&readid_file)
                .stdin(Stdio::piped())
                .stdout(Stdio::piped())
                .spawn()?;

            let fastq_path = fastq.clone();
            let stdin = child.stdin.take().unwrap();
            let writer = thread::spawn(move || -> io::Result<()> {
                let mut decoder = GzDecoder::new(File::open(&fastq_path)?);
                let mut w = io::BufWriter::new(stdin);
                io::copy(&mut decoder, &mut w)?;
                Ok(())
            });

            io::copy(&mut child.stdout.take().unwrap(), outfile)?;
            writer.join().expect("gz writer thread panicked")?;
            child.wait()?;
        } else {
            let mut child = Command::new("seqtk")
                .arg("subseq")
                .arg(fastq)
                .arg(&readid_file)
                .stdout(Stdio::piped())
                .spawn()?;
            io::copy(&mut child.stdout.take().unwrap(), outfile)?;
        }
        Ok(())
    };

    process_seqtk(&fastq_files[0], &mut outfile1)?;
    // Process reverse reads
    if let Some(outfile2) = &mut outfile2 {
        process_seqtk(&fastq_files[1], outfile2)?;
    }

    Ok(())
}
