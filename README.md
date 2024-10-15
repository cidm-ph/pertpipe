# pertpipe
Pipeline for full analysis of Bordetella pertussis, particularly for the detection of macrolide resistance &amp; vaccine antigen mutations.

## Development
The initial development of pertpipe was designed for culture-based sequencing methods, ensuring accurate detection of 23S rRNA mutations and the vaccine antigens. Pertpipe primarily reports information about the ptxP, ptxA-E, prn, fim2, fim3 types and any mutations occuring the V domain of the 23S rRNA. As CIDM-PH transition to the use of the hybridisation capture-based method for _B. pertussis_, the pipeline has been adapted to take metagenomic options.

## Validation
The pipeline was validated using 30 paired sequences.

# Installation
PIP (Recommended)
```
pip install git+https://github.com/cidm-ph/pertpipe
```
GITHUB
```
git clone https://github.com/cidm-ph/pertpipe.git
```
Then you can use setup.py
```
python setup.py
```

## Usage
Example usage

Reads only
```
pertpipe --R1 $PATH/$R1.fq.gz --R2 $PATH/$R2.fq.gz --outdir $OUTDIR
```

Fasta only
```
pertpipe --fasta $PATH/example.fasta --outdir $OUTDIR
```

Metagenomics
Fasta only
```
pertpipe --R1 $PATH/$R1.fq.gz --R2 $PATH/$R2.fq.gz --outdir $OUTDIR --meta
```

Its also possible to skip assembly, annotation or some mapping steps by supplying the entire output folder of Spades, Megahit or Prokka as a folder in the output directory.
Prokka is also not a necessary step to perform for the analysis but is useful in general for drawing trees etc.

FLAGS

```
--outdir, -o [PATH]             optional folder to write output files to
--R1 [PATH]                     R1 fastq of sample (can be gzipped files)
--R2 [PATH]                     R2 fastq of sample (can be gzipped files)
--fasta, -f [PATH]              optional fasta file which will skip spades, can be used in conjunction with --longread.             
--meta, -m                      optional flag to trigger metagenomics mode, which utilises meta flags of prokka and metaSPades
--version, -v                   print version
```

## Output
- A SPades, metaSPades or megahit folder (spades/megahit)
- A Prokka folder (prokka)
- An analysis folder containing several intermediate files (analysis)
- "vir_res.tsv", which is the file containing the final result.
- A log file which you can use to view intermediate decisions.

## Limitations
- If using Metagenomics, it is important to remember that the 23S rRNA is quite a conservative region across many species, therefore commensals can play a part in creating "noise" and as a result, induce false SNP mutations in this region. The pipeline will do its best to remove as much of this noise as possible but it cannot guarantee that SNPs detected in the 23S rRNA confer resistance to macrolides (except A2037G).

## Dependencies
- SPades
- Prokka
- Samtools
- minimap2
- megahit
- Kallisto


