# telo_tools

```
python -m venv telo_tools
source telo_tools/bin/activate

git clone https://github.com/touala/telo_tools
pip install -r telo_tools/requirements.txt
```

```
python telo_tools/evaluate_telomere_assembly.py -i <input.bam> -t <telofinder_annotation> -o <output.tsv> -m <int_margin>
python telo_tools/extract_telomere_sequence.py -i <input.bam> -t <telofinder_annotation> -o <output.tsv> -c -s
```
