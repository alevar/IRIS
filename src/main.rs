extern crate clap;
extern crate rust_htslib;

use std::fs;
use std::io::Write;
use std::fs::File;
use std::io::{BufRead, BufReader};

use clap::{Command,Arg, ArgAction};
use rust_htslib::bam;
use rust_htslib::bam::Read;

fn extract_genes(fname: &str) {

    // Open the GTF file
    if let Ok(file) = File::open(fname) {
        let reader = BufReader::new(file);

        for line in reader.lines() {
            if let Ok(line) = line {
                if line.starts_with('#') {
                    continue; // Skip comment lines
                }

                let lcs: Vec<&str> = line.split('\t').collect();
                if lcs.len() < 9 {
                    continue; // Skip lines with insufficient fields
                }

                let feature_type = lcs[2];
                if feature_type != "exon" && feature_type != "transcript" {
                    continue; // Skip non-exon and non-transcript lines
                }

                let attributes: Vec<&str> = lcs[8].split(';').collect();
                let mut gene_id = "";
                for attr in attributes {
                    let parts: Vec<&str> = attr.trim().split(' ').collect();
                    if parts.len() >= 2 && parts[0] == "gene_id" {
                        gene_id = parts[1].trim_matches('"');
                        break;
                    }
                }

                let seqid = lcs[0].to_string();
                let strand = lcs[6].chars().next().unwrap_or('+');
                let start = lcs[3].parse::<usize>().unwrap_or(0);
                let end = lcs[4].parse::<usize>().unwrap_or(0) + 1; // Exclusive end

                println!("{} {} {} {} {}", seqid, start, end, strand, gene_id);
            }
        }
    }
}

// the entry point for the software
fn run(i1: &String, i2: &String, a1: &String, a2: &String, o: &String) {
    println!("Input file 1: {}", i1);
    extract_genes(a1);
}

fn main() {
    let matches = Command::new("IRIS")
        .version("0.1.0")
        .author("Ales Varabyou")
        .about("Inference of RetroViral Integration Sites")
        .arg(
            Arg::new("i1")
            .long("i1")
            .help("File containing mapping of reads to genome #1.")
        )
        .arg(
            Arg::new("i2")
            .long("i2")
            .help("File containing mapping of reads to genome #2.")
        )
        .arg(
            Arg::new("a1")
            .long("a1")
            .required(false)
            .help("GTF file containing gene annotations for genome #1.")
        )
        .arg(
            Arg::new("a2")
            .long("a2")
            .required(false)
            .help("GTF file containing gene annotations for genome #2.")
        )
        .arg(
            Arg::new("output")
            .short('o')
            .long("output")
            .help("Output file name."),

        )
        .after_help("--help or -h")
        .get_matches();

    // Extract and use the value of the "input" argument
    let input1: &String = matches.get_one("i1").unwrap();
    let input2: &String = matches.get_one("i2").unwrap();
    let gtf1: &String = matches.get_one("a1").unwrap();
    let gtf2: &String = matches.get_one("a2").unwrap();
    let output: &String = matches.get_one("output").unwrap();

    run(input1, input2, gtf1, gtf2, output);
}
