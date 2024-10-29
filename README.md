# telo_tools

Accessories scripts for `Nanotiming` manuscript [Theulot et al., 2024](https://doi.org/10.1101/2024.07.05.602252). See the main `Nanotiming` repository [here](https://github.com/LacroixLaurent/NanoTiming).

```
#######
# Installation
#######

python -m venv telo_tools
source telo_tools/bin/activate

git clone https://github.com/touala/telo_tools
pip install -r telo_tools/requirements.txt
```

```
#######
# Default usage
#######

# Compute basic metrics to evaluate the quality of telomere in an assembly
#   First, report alignment and soft-clipped lengths for each read alignemnts (including non-primary)
#   Second, compute basic statistics for primary alignments only
python telo_tools/evaluate_telomere_assembly.py -i <input.bam> -t <telofinder_annotation> -o <output.tsv> -m <int_margin>
cat <output.tsv> | awk '{if($3 == 0 || $3 == 16){nb_reads+=1; size_aln+=$5; size_sc+=$6+$7; }} END{print(nb_reads, size_aln, size_aln/nb_reads, size_sc, size_sc/nb_reads)}'

# Extract all telomeric sequences for each chromsome end, see main Nanotiming repository for further filtering
python telo_tools/extract_telomere_sequence.py -i <input.bam> -t <telofinder_annotation> -o <output.tsv> -c -s
```
