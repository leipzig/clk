```
snakemake --no-shared-fs --default-remote-provider S3 --default-remote-prefix clk  \
clk/SRP091981/untreated_vs_0.05.manifest.txt clk/SRP091981/untreated_vs_0.1.manifest.txt clk/SRP091981/untreated_vs_0.5.manifest.txt clk/SRP091981/untreated_vs_1.0.manifest.txt
```

Here I reproduced the differential isoform usage induced by T3 using rMATS, rMATS-ISO, SUPPA2, and (the original paper used GSNAP/MISO/Cufflinks)
<img width="689" alt="Screen Shot 2021-04-21 at 12 47 24 PM" src="https://user-images.githubusercontent.com/147991/115591081-f8e8e380-a28e-11eb-85f3-fa8d4c89b5bd.png">

|                          | Original implementation | YX        | EL     | AN          |
|--------------------------|-------------------------|-----------|--------|-------------|
| Aligner                  | GSNAP                   | STAR      | Salmon |             |
| Genome                   | GRCh37                  | GRCh38    | GRCh38 | GRCh38      |
| Differential Splicing    | MISO                    | rMATS     | SUPPA2 |             |
| Isoform Detection        |                         | rMATS-ISO |        | MINTIE      |
| Conjoined gene detection | deFuse                  |           |        | STAR-FUSION |

Reviewers are free to alter whatever they choose in order to assess the validity of this study. 

Some suggestions include, but are not limited to...
