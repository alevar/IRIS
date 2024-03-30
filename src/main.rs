extern crate clap;
use clap::{Command, Arg, ArgAction};

mod iris {
    use std::collections::{HashMap, HashSet};
    use std::fs;
    use std::hash::Hash;
    use std::io::{self, BufRead, Write};
    use bio::utils::Interval;
    use bio::data_structures::interval_tree::ArrayBackedIntervalTree;

    fn extract_attributes(attribute_str: &str) -> std::collections::HashMap<String, String> {
        let mut attrs_dict = std::collections::HashMap::new();
        let attrs: Vec<&str> = attribute_str.trim().trim_end_matches(';').split(';').map(|s| s.trim()).collect();
    
        for attr in attrs {
            let parts: Vec<&str> = attr.trim_matches('"').splitn(2, " \"").collect();
            if parts.len() == 2 {
                let key = parts[0].to_string();
                let value = parts[1].trim_matches('"').to_string();
                attrs_dict.insert(key, value);
            }
        }
    
        attrs_dict
    }

    pub fn extract_genes(fname: &str) -> Result<(HashMap<(String, char), ArrayBackedIntervalTree<u32, String>>, HashMap<(String, char), ArrayBackedIntervalTree<u32, String>>), io::Error> {
        let mut exon_map: HashMap<(String, char), ArrayBackedIntervalTree<u32, String>> = HashMap::new();
        let mut intron_map: HashMap<(String, char), ArrayBackedIntervalTree<u32, String>> = HashMap::new();
    
        let file = fs::File::open(fname)?;
        let reader = io::BufReader::new(file);
    
        let mut gene_id = String::new();
        let mut gene_name = String::new();
        let mut transcript_id = String::new();
    
        let mut prev_exon: Option<(u32, u32)> = None;
    
        for line in reader.lines() {
            let line = line?;
            if line.starts_with('#') {
                continue; // Skip comment lines
            }
    
            let lcs: Vec<&str> = line.split('\t').collect();
            if lcs.len() < 9 {
                continue; // Skip lines with insufficient fields
            }
    
            let feature_type = lcs[2];
            if feature_type != "exon" && feature_type != "transcript" {
                continue; // Skip non-exon lines
            }
    
            let attributes: HashMap<String, String> = extract_attributes(lcs[8]);
    
            let tid = match attributes.get("transcript_id") {
                Some(value) => value.clone(),
                None => return Err(io::Error::new(io::ErrorKind::Other, "All entries are expected to have transcript_id attribute")),
            };
    
            match feature_type {
                "transcript" => {
                    transcript_id = tid.clone();
                    prev_exon = None;
                    gene_id = match attributes.get("gene_id") {
                        Some(value) => value.to_string(),
                        None => return Err(io::Error::new(io::ErrorKind::Other, "All transcript entries are expected to have gene_id attribute")),
                    };
                    gene_name = match attributes.get("gene_name") {
                        Some(value) => value.to_string(),
                        None => return Err(io::Error::new(io::ErrorKind::Other, "All transcript entries are expected to have gene_name attribute")),
                    };
                }
                "exon" => {
                    if tid != transcript_id {
                        return Err(io::Error::new(io::ErrorKind::Other, format!("Exon entry with transcript_id {} does not match the last seen transcript_id {}. Make sure the GTF file is correctly formatted and sorted. You can use gffread -T to standardize your annotation", tid, transcript_id)));
                    }
    
                    let seqid = lcs[0].to_string();
                    let strand = lcs[6].chars().next().unwrap_or('+');
                    let start: u32 = lcs[3].parse().unwrap_or(0);
                    let end: u32 = lcs[4].parse().unwrap_or(0) + 1; // Exclusive end
    
                    let exon_entry = exon_map.entry((seqid.clone(), strand)).or_insert_with(|| ArrayBackedIntervalTree::new());
                    exon_entry.insert(Interval::new(start..end).unwrap(), gene_id.clone());
    
                    if let Some((prev_start, prev_end)) = prev_exon {
                        let intron_entry = intron_map.entry((seqid.clone(), strand)).or_insert_with(|| ArrayBackedIntervalTree::new());
                        intron_entry.insert(Interval::new(prev_end..start).unwrap(), gene_id.clone());
                    }
    
                    prev_exon = Some((start, end));
                }
                _ => continue,
            }
        }
    
        // Perform indexing of the interval trees
        for tree in exon_map.values_mut() {
            tree.index();
        }
        for tree in intron_map.values_mut() {
            tree.index();
        }
    
        Ok((exon_map, intron_map))
    }

    pub fn extract_donor_acceptor(intron_map: &HashMap<(String,char), ArrayBackedIntervalTree<u32, String>>) -> Result<(HashMap<String, HashMap<u32,(String,char)>>, HashMap<String, HashMap<u32,(String,char)>>),io::Error> {
        let mut donors: HashMap<String, HashMap<u32,(String,char)>> = HashMap::new();
        let mut acceptors: HashMap<String, HashMap<u32,(String,char)>> = HashMap::new();

        for ((seqid,strand), intron_tree) in intron_map.iter() {
            donors.insert(seqid.clone(), HashMap::new());
            acceptors.insert(seqid.clone(), HashMap::new());

            for (intron) in intron_tree {
                let start = intron.interval().start;
                let end = intron.interval().end;

                donors.get_mut(seqid).unwrap().insert(start, (intron.data().clone(), strand.clone()));
                acceptors.get_mut(seqid).unwrap().insert(end, (intron.data().clone(), strand.clone()));
            }
        }

        Ok((donors, acceptors))
    }
}

// the entry point for the software
fn run(i1: &String, i2: &String, a1: &String, a2: &String, o: &String) {
    let (a1_exon_map, a1_intron_map) = match iris::extract_genes(a1) {
        Ok(res) => res,
        Err(err) => {
            eprintln!("Error: {}", err);
            std::process::exit(1);
        }
    };

    let (a2_exon_map, a2_intron_map) = match iris::extract_genes(a2) {
        Ok(res) => res,
        Err(err) => {
            eprintln!("Error: {}", err);
            std::process::exit(1);
        }
    };

    let (a1_donors, a1_acceptors) = match iris::extract_donor_acceptor(&a1_intron_map) {
        Ok(res) => res,
        Err(err) => {
            eprintln!("Error: {}", err);
            std::process::exit(1);
        }
    };

    let (a2_donors, a2_acceptors) = match iris::extract_donor_acceptor(&a2_intron_map) {
        Ok(res) => res,
        Err(err) => {
            eprintln!("Error: {}", err);
            std::process::exit(1);
        }
    };
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
