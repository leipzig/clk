```
wget https://s3.amazonaws.com/panorama-public/manifest.txt && wget --input-file manifest.txt --base https://s3.amazonaws.com/panorama-public/
```

or preferably...
```
aws s3 sync s3://panorama-public/ .
```

to run
```
time docker run -v $PWD:/rMATS-ISO-master/work quay.io/leipzig/rmats-iso python ./rMATS-ISO.py --in-gtf work/level_1_protein_coding_genes.gtf --in-bam work/untreatedvslowdose.manifest.txt -o /work
```

above takes about 7 minutes with `level_1_protein_coding_genes.gtf` which has 5573 genes and 3456 transcripts, 33158 exons



with 
```
time docker run -v $PWD:/rMATS-ISO-master/work quay.io/leipzig/rmats-iso python ./rMATS-ISO.py --in-gtf work/Homo_sapiens.GRCh38.87.gtf --in-bam work/untreatedvslowdose.manifest.txt -o /work
```


gtf sizes
```
for f in *gtf; do echo $f; cut -f3 $f | grep -v '#' | sort | uniq -c; done 
gencode.v28.basic.annotation.gtf
 555199 CDS
 701383 exon
  58381 gene
     97 Selenocysteine
  56753 start_codon
  56641 stop_codon
 100985 transcript
 155098 UTR
genes.gtf
 864338 CDS
1031317 exon
Homo_sapiens.GRCh38.87.gtf
 703935 CDS
1182163 exon
 142387 five_prime_utr
  58051 gene
    117 Selenocysteine
  82657 start_codon
  74244 stop_codon
 133938 three_prime_utr
 198002 transcript
level_1_genes.gtf
  30336 CDS
  50534 exon
  16725 gene
      6 Selenocysteine
   3457 start_codon
   3465 stop_codon
  14524 transcript
   9533 UTR
level_1_protein_coding_genes.gtf
  30322 CDS
  33158 exon
   5573 gene
      6 Selenocysteine
   3455 start_codon
   3464 stop_codon
   3456 transcript
   9529 UTR
```