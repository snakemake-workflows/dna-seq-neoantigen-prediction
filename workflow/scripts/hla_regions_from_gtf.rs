//! This is a regular crate doc comment, but it also contains a partial
//! Cargo manifest.  Note the use of a *fenced* code block, and the
//! `cargo` "language".
//!
//! ```cargo
//! cargo-features = ["edition2021"]
//! [dependencies]
//! bio = { version = "0.41.0" }
//! ```
use bio::io::{gff, bed};

use std::fs::File;
use std::path::PathBuf;
use std::collections::HashSet;
use std::io::{BufRead, BufReader};
use std::error::Error;


fn main() -> Result<(), Box<dyn Error>> {

    snakemake.redirect_stderr(&snakemake.log[0])?;

    let alleles_file = BufReader::new(File::open(PathBuf::from(&snakemake.input.allele_names))?);
    let allele_names: HashSet<String> = alleles_file.lines().map(|line| line.unwrap()).collect();

    let mut gtf_reader = gff::Reader::from_file(PathBuf::from(&snakemake.input.gtf), gff::GffType::GTF2)?;

    let mut bed_writer = bed::Writer::to_file(PathBuf::from(&snakemake.output[0]))?;

    for r in gtf_reader.records() {
        let record = r?;
        if record.feature_type() == "gene" {
            let attr = record.attributes();
            if let Some(name) = attr.get("gene_name") {
                if allele_names.contains(name) {
                    let mut bed_record = bed::Record::new();
                    bed_record.set_chrom(record.seqname());
                    bed_record.set_start(*record.start());
                    bed_record.set_end(*record.end());
                    bed_record.set_name(name);
                    // write out bed record
                    bed_writer.write(&bed_record)?;
                }
            } 
        }
    }

    Ok(())
}