use std::io;
use std::ops::Range;
use std::collections::{HashMap, HashSet};

// find_genomes needs these
extern crate csv;
#[macro_use] extern crate log;
extern crate env_logger;
extern crate serde;
#[macro_use] extern crate serde_derive;

// Import find_genomes
mod find_genomes;
use find_genomes::{find_genomes, Gene, Genome};

fn main() {
    // Initialize the logger.
    // To see warnings, use `RUST_LOG=warn /path/to/intersections`
    // The `.ok()` specifies that the program can continue, even if this fails.
    env_logger::init().ok();
    // The data that will appear in the output file.
    // Maps a gene to its per genome counts.
    // `HashMap<a, b>` is a map/dictionary: given `a`, get `b`.
    // `Vec<u32>` is a list of positive integers (the integers are 32 bits each).
    let mut data: HashMap<Gene, Vec<u32>> = HashMap::new();
    // We collect the genomes here so we can see how many there are.
    // Normally, iterators are generated on the fly, so the total count is not known.
    // The "." is the directory to collect genomes from.
    // `expect` specifies how to handle an IO error.
    // Upon an error, `expect` will print an error than exit (known as panicking).
    // `expect` is like `unwrap`, except we add more information to the error message.
    // A Result is a value that may be an error.
    // An `io::Result<Genome>` may be either `Ok(Genome)` or `Err(io::Error)`.
    // We'll call `expect` on it later to turn it into just a `Genome`.
    // Just like the `expect` call here, it will panic (crash) if it's an error.
    // For more information, see:
    // https://doc.rust-lang.org/std/result/enum.Result.html#method.expect
    // https://doc.rust-lang.org/std/result/
    // https://rustbyexample.com/error.html
    let genomes: Vec<io::Result<Genome>> = find_genomes(".").expect("Failed to list genomes").collect();
    let genome_count = genomes.len();
    // We'll fill this list up later
    // It's faster to fill it up later because we're already looping through it.
    // We could just do `Vec::new()`, the `with_capacity` is a slight optimization.
    let mut genome_names: Vec<String> = Vec::with_capacity(genome_count);
    // Iterate over the genomes.
    // Using the `enumerate` method also tells us the index of the genome (`i`).
    // We'll need it later on.
    for (i, genome) in genomes.into_iter().enumerate() {
        // The `expect` method specifies how to handle an IO error
        // See the previous `expect` for find_genomes earlier for more details
        let genome = genome.expect("Failed to read genomes");
        // Add this genome name to the list of names
        genome_names.push(genome.name);
        // Create a map of accn to a list of sequences
        let mut sequences: HashMap<Vec<u8>, Vec<Range<usize>>> = HashMap::new();
        for blast_line in genome.blast_iter {
            // This line does two things. First, it specifies how to handle an error.
            // It also splits up the blast line into the accn and sequence.
            let (accn, range) = blast_line.expect("Failed to read BLAST file");
            // TODO add to the sequences map
	    let sequence = sequences.entry(accn).or_insert_with(Vec::new);
            if range.end > sequence_counts.len() {
                sequences.resize(range.end, 0);
            }
	    for index in range {
                sequences[index] += 1;
            }
            // hint: later on we use the `entry` API - could that be useful here?
            // hint: the default should be `Vec::new()`, meaning no sequences are present
            // https://doc.rust-lang.org/std/collections/struct.HashMap.html#method.entry
        }
        // A `HashSet` is basically an unordered list
        // This one will keep track of missing accns, as to avoid spamming the logs with them
        let mut accns_missing: HashSet<Vec<u8>> = HashSet::new();
        for gff_line in genome.gff_iter {
            // See previous `blast_line.expect` line for explanation
            let (accn, gene, range) = gff_line.expect("Failed to read GFF file");
            let mut count = 0;
            // The `if let` structure checks if the right hand of the `=` matches the left side.
            // The `&` indicated that we are only temporarily using the gene's accn.
            if let Some(sequences) = sequences.get(&accn) {
                // This is executed if it does match, and there are sequences with this accn.
                for sequence in sequences {
                    // TODO check if we should add to the count
                    
                    // Hint: `sequence` is of type range
                    // Hint: You'll also want to use `range`, the range of this gene
                    // Hint: ranges have start and end fields
                    // https://doc.rust-lang.org/std/ops/struct.Range.html
                }
            } else {
                // If this is the first time we've encountered this accn, log it.
                // `insert` will return true if the value isn't already in the set.
                // https://doc.rust-lang.org/std/collections/struct.HashSet.html#method.insert
                // We need to clone the accn, as we'll be using it later and `insert` consumes it.
                // One copy of the accn gets put in the HashSet and the other is logged.
                // There's a faster way to do this without clone, but it's a lot more complicated.
                // Speed isn't important here as this should theoretically never happen.
                if accns_missing.insert(accn.clone()) {
                    // For logging, attempt to turn the accn into a human readable string
                    // https://doc.rust-lang.org/std/string/struct.String.html#method.from_utf8_lossy
                    let accn = String::from_utf8_lossy(&accn);
                    // The {} will be substituted in for the first argument, accn
                    warn!("Encountered an accn with a gene but no sequence: {}", accn);
                }
            }
            // Get the entry for the gene in the data.
            // If it doesn't exist, put in an array of 0s.
            // We use `entry` here, and not `data[gene]`, because it might not exist yet.
            // Doing `data[gene]` would crash if data didn't have a key `gene`.
            // The `||` defines a "closure", basically a mini-function.
            // This function always returns a vector of length `genome_count` filled with 0s.
            // The `vec!` is short hand for that.
            // We could do `.or_insert(vec![0; genome_count])` instead, but that'd be slower.
            // If we did that, we'd have to create the `Vec` each time, even if we didn't use it.
            let gene_counts = data.entry(gene).or_insert_with(|| vec![0; genome_count]);
            // TODO add the count to correct field in gene_counts
            // Hint: `gene_counts` is of type `Vec<u32>` (`u32` is a positive integer)
            // Hint: `i` is the index of the genome, as defined up above
        }
    }
    // Start by writing the header
    print!("Name\tProduct");
    for name in genome_names {
        // The {} is substituted in for the next argument passed to `print!`.
        // \t is a tab.
        print!("\t{}", name);
    }
    // Begin the next line.
    println!("");
    // We've collected the data, now loop through it
    for (gene, genome_counts) in data {
        // This statement is a bit fancy
        // It checks if the gene specifies a name
        if let Some(name) = gene.name {
            // This is executed if it specifies a name
            print!("{}\t", name);
        } else {
            // This is executed if doesn't specify a name
            print!("Unknown\t");
        }
        // Genes always have a product
        print!("{}", gene.product);
        // TODO print out the data
        // Hint: `genome_counts` is an array of positive integers
        // Hint: It lines up with the order of the genomes in the headers
        
        // Begin the next line
        println!("");
    }
}
