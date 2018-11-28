To run standard comparisons:
```
snakemake --no-shared-fs --default-remote-provider S3 --default-remote-prefix panorama-clk-repro  \
panorama-clk-repro/SRP091981/untreated_vs_0.05.manifest.txt panorama-clk-repro/SRP091981/untreated_vs_0.1.manifest.txt panorama-clk-repro/SRP091981/untreated_vs_0.5.manifest.txt panorama-clk-repro/SRP091981/untreated_vs_1.0.manifest.txt
```