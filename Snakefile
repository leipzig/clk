from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
from snakemake.utils import R
from utils import metautils
import yaml
import boto3
import re

configfile: "config.yaml"

PROJECT_BUCKET = 'clk-splicing/SRP091981'


LOCAL_SCRATCH = "/scratch"
RAWDIR="SRP091981"
PROCESSDIR="process"
SRAFILES = [line.rstrip() for line in open("metadata/SraAccList.txt")]
ILLUMINA_SRA = metautils.illuminaRuns()
PACBIO_SRA = metautils.pacbioRuns()
TMPDIR = os.getcwd()

rule onemini:
    input: RAWDIR+"/SRR5009429_sam_novel.gtf", RAWDIR+"/SRR5009429.lr2rmats.log", RAWDIR+"/SRR5009429.filtered.bam"

rule justone:
    input: RAWDIR+"/SRR5009515.Aligned.sortedByCoord.out.bam"

rule bigboy:
    input: RAWDIR+"/SRR5009517.Aligned.sortedByCoord.out.bam"

   
rule s3_illumina_files:
    input: expand(RAWDIR+"/{sampleids}_{pair}.fastq.gz", sampleids=ILLUMINA_SRA, pair=[1,2])

rule illumina_align:
    input: expand(RAWDIR+"/{sampleids}.Aligned.sortedByCoord.out.bam", sampleids=ILLUMINA_SRA)

rule illumina_index:
    input: expand(RAWDIR+"/{sampleids}.Aligned.sortedByCoord.out.bam.bai", sampleids=ILLUMINA_SRA)


rule s3_pacbio_files:
    input: expand(RAWDIR+"/{sampleids}.fastq.gz", sampleids=PACBIO_SRA)

rule pacbio_align:
    input: expand(RAWDIR+"/{sampleids}.sam", sampleids=PACBIO_SRA)

rule pacbio_lr2rmats:
    input: expand(RAWDIR+"/{sampleids}_sam_novel.gtf", sampleids=PACBIO_SRA)

rule allfiles:
    input: expand(RAWDIR+"/{sampleids}_{pair}.fastq.gz", sampleids=ILLUMINA_SRA, pair=[1,2]), expand(RAWDIR+"/{sampleids}.fastq.gz", sampleids=PACBIO_SRA)

rule fetchpair_from_aws:
    output: pair1 = RAWDIR+"/{accession}_1.fastq.gz",
            pair2 = RAWDIR+"/{accession}_2.fastq.gz"
    run:
        s3pair1 = metautils.st.loc[metautils.st['Run'] == wildcards.accession]['Pair1Filename'].to_string(index=False).replace(' ','')
        s3pair2 = metautils.st.loc[metautils.st['Run'] == wildcards.accession]['Pair2Filename'].to_string(index=False).replace(' ','')
        print(s3pair1)
        if wildcards.accession in metautils.pacbioRuns():
            #pair1 correspond to the bas.h5
            shell("wget https://clk-splicing.s3.amazonaws.com/SRP091981/{0}/{1}".format(wildcards.accession,s3pair1))
            shell("bash5tools.py {0} --outFilePrefix {1}".format(s3pair1,wildcards.accession))
        	#shell("bedtools bamtofastq -i {0} -fq {1} -fq2 {2}".format(output.pair1,output.pair2))
        else:
            shell("wget -O {2} https://clk-splicing.s3.amazonaws.com/SRP091981/{0}/{1}".format(wildcards.accession,s3pair1,output.pair1))
            shell("wget -O {2} https://clk-splicing.s3.amazonaws.com/SRP091981/{0}/{1}".format(wildcards.accession,s3pair2,output.pair2))
       

# Let's not do this - I just did a bulk transfer of SRP091981
# rule fetchpair_from_sra:
#     output: RAWDIR+"/{accession}_1.fastq.gz", RAWDIR+"/{accession}_2.fastq.gz"
#     run:
#         print("fastq-dump --split-3 --gzip {0} -O {1}".format(wildcards.accession,RAWDIR))
#         shell("fastq-dump --split-3 --gzip {0} -O {1}".format(wildcards.accession,RAWDIR))
#         #shell("fasterq-dump --split-3 {0} -O {1}".format(wildcards.accession,RAWDIR))
#         #shell("gzip -1 {0}/{1}_1.fastq {0}/{1}_2.fastq".format(RAWDIR,wildcards.accession))

rule getasingleton:
    output: RAWDIR+"/{accessiondf -h }.fastq.gz"
    run:
        print("fastq-dump --split-3 --gzip {0} -O {1}".format(wildcards.accession,RAWDIR))
        shell("fastq-dump --split-3 --gzip {0} -O {1}".format(wildcards.accession,RAWDIR))

#this is a more casual metadata file than SRA provides through pysradb
rule metadata:
    output: "metadata.csv"
    shell:
        """
        wget -O metadata.csv 'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=SRP091981'
        """

# rule gzip:
#     input: "{file}"
#     output: "{file}.gz"
#     shell: "gzip {input}"

#we lock this down to gz for disambiguate
rule upload:
    input: "{file}.gz"
    output: PROJECT_BUCKET+"/{file}.gz"
    run:
          shell("mv {input} {output} && rm -f {input}")

rule getStarRefs:
    output:
        "GRCh38_star/Genome",
        "GRCh38_star/genomeParameters"
    shell:
        """
        aws s3 sync s3://clk-splicing/refs/GRCh38/Sequence/STARIndex/ GRCh38_star/
        """

rule star_align:
    input: RAWDIR+"/{sample}_1.fastq.gz", RAWDIR+"/{sample}_2.fastq.gz"
    output: RAWDIR+"/{sample}.Aligned.sortedByCoord.out.bam",
            RAWDIR+"/{sample}.Log.final.out",
            RAWDIR+"/{sample}.Log.out",
            RAWDIR+"/{sample}.Log.progress.out",
            RAWDIR+"/{sample}.SJ.out.tab"
    threads: 8
    params: bytes = lambda wildcards: metautils.getECS(wildcards.sample,'bytes','STAR'),
            mb = lambda wildcards: metautils.getECS(wildcards.sample,'mb','STAR')
    shell: """
            STAR --runMode alignReads \
                 --outSAMtype BAM SortedByCoordinate \
                 --limitBAMsortRAM {params.bytes} \
                 --readFilesCommand zcat \
                 --outFilterType BySJout   --outFilterMultimapNmax 20 \
                 --outFilterMismatchNmax 999   --alignIntronMin 25  \
                 --alignIntronMax 1000000   --alignMatesGapMax 1000000 \
                 --alignSJoverhangMin 8   --alignSJDBoverhangMin 5 \
                 --sjdbGTFfile GRCh38_star/genes.gtf \
                 --genomeDir GRCh38_star \
                 --runThreadN {threads} \
                 --outFileNamePrefix SRP091981/{wildcards.sample}.  \
                 --readFilesIn  SRP091981/{wildcards.sample}_1.fastq.gz SRP091981/{wildcards.sample}_2.fastq.gz
            
            samtools index SRP091981/{wildcards.sample}.Aligned.sortedByCoord.out.bam
            """

# minimap mapping for long reads
rule minimap_map:
    input:
        "SRP091981/{sample}.fasta"
    output:
        "SRP091981/{sample}.sam"
    threads:
        config["minimap_map"]["threads"]
    log:
        "logs/minimap_map/{sample}.log"
    benchmark:
        "benchmark/{sample}.minimap.benchmark.txt"
    params:
        minimap=config["exe_files"]["minimap2"],
        mb = lambda wildcards: metautils.getECS(wildcards.sample,'mb','minimap')
    shell:
        """
        export sample="{wildcards.sample}""
        export project="SRP091981"
        sh scripts/minimap.sh
        """

rule generate_two_way_manifest:
    output: manifest=RAWDIR+"/{sample1,[a-z0-9.-]+}_vs_{sample2,[a-z0-9.-]+}.manifest.txt"
    run:
        metautils.twoSampleComparisonManifest(wildcards.sample1,wildcards.sample2,output.manifest)



rule run_rmatsiso_from_bam:
    input: bam=RAWDIR+"/{sample}.Aligned{ext}.bam", bai=RAWDIR+"/{sample}.Aligned{ext}.bam.bai"
    output: RAWDIR+"/{sample}.Aligned{ext}.bam.IsoExon", RAWDIR+"/{sample}.Aligned{ext}.bam.IsoMatrix"
    params: bytes = lambda wildcards: metautils.getECS(wildcards.sample,'bytes','IsoModule'),
            mb = lambda wildcards: metautils.getECS(wildcards.sample,'mb','IsoModule'),
            gtf = "level_1_protein_coding_genes.gtf"
    shell: """
            export sample="{wildcards.sample}.Aligned{wildcards.ext}"
            export project="SRP091981"
            export bytes="{params.bytes}"
            export gtf="{params.gtf}"
            sh scripts/isomodule.sh
            """

rule run_rmatsiso_from_manifest:
    input: untreated=lambda wildcards: metautils.getBamsFromSampleName(metautils.getfulldosagename(wildcards.sample1),include_s3=RAWDIR),
           treated=lambda wildcards: metautils.getBamsFromSampleName(metautils.getfulldosagename(wildcards.sample2),include_s3=RAWDIR),
           manifest=RAWDIR+"/{sample1}_vs_{sample2}.manifest.txt"
    output: manifest=RAWDIR+"/{sample1}_vs_{sample2}/done"
    params: bytes = lambda wildcards: metautils.getECS('foo','bytes','IsoModule'),
            mb = lambda wildcards: metautils.getECS('foo','mb','IsoModule'),
            gtf = "gencode.v28.annotation.gtf",
            jobname = lambda wildcards: re.sub('\.','',wildcards.sample1+'_'+wildcards.sample2)
    shell: """
            export untreated="{input.untreated}"
            export treated="{input.treated}"
            export manifest="{input.manifest}"
            export comparison="{wildcards.sample1}_vs_{wildcards.sample2}"
            export project="SRP091981"
            export bytes="{params.bytes}"
            export gtf="{params.gtf}" 
            sh scripts/isomanifest.sh
            """

rule run_rmatsturbo_from_manifest:
    input: untreated=lambda wildcards: metautils.getBamsFromSampleName(wildcards.sample1,include_s3=RAWDIR),
           treated=lambda wildcards: metautils.getBamsFromSampleName(wildcards.sample2,include_s3=RAWDIR),
           manifest=RAWDIR+"/{sample1}_vs_{sample2}.manifest.txt"
    output: manifest=RAWDIR+"-turbo/{sample1,[a-z0-9.-]+}_vs_{sample2,[a-z0-9.-]+}/done"
    params: bytes = lambda wildcards: metautils.getECS('foo','bytes','IsoModule'),
            mb = lambda wildcards: metautils.getECS('foo','mb','IsoModule'),
            gtf = "gencode.v28.annotation.gtf",
            reftx = "GRCh38_star",
            jobname = lambda wildcards: re.sub('\.','',wildcards.sample1+'_'+wildcards.sample2)
    shell: """
            export untreated="{input.untreated}"
            export treated="{input.treated}"
            export manifest="{input.manifest}"
            export comparison="{wildcards.sample1}_vs_{wildcards.sample2}"
            export project="SRP091981-turbo"
            export bytes="{params.bytes}"
            export gtf="{params.gtf}" 
            export reftx="{params.reftx}"
            export gtf="{params.gtf}"
            sh scripts/turbomanifest.sh
            """

rule all_sashimi:
    input: expand(RAWDIR+"-sashimi/{sample1}_vs_{sample2}/done",sample1="untreated",sample2=['0.05','0.1','0.5','1.0','treated'])

#snakemake --force --no-shared-fs --default-remote-provider S3 --default-remote-prefix panorama-clk-repro  --cluster runners/slurm_scheduler.py --cluster-status runners/eric_status.py  -j 50 --cluster-config runners/slurm_cluster_spec.yaml panorama-clk-repro/SRP091981-sashimi/untreated_vs_treated/done panorama-clk-repro/SRP091981-sashimi/untreated_vs_0.1/done  panorama-clk-repro/SRP091981-sashimi/untreated_vs_0.5/done panorama-clk-repro/SRP091981-sashimi/untreated_vs_1.0/done 
rule run_rmatssashimi_from_manifest:
    input: untreated=lambda wildcards: metautils.getBamsFromSampleName(wildcards.sample1,include_s3=RAWDIR),
           treated=lambda wildcards: metautils.getBamsFromSampleName(wildcards.sample2,include_s3=RAWDIR),
           manifest=RAWDIR+"/{sample1}_vs_{sample2}.manifest.txt",
           rmats=RAWDIR+"-turbo/{sample1}_vs_{sample2}/done"
    output: manifest=RAWDIR+"-sashimi/{sample1,[a-z0-9.-]+}_vs_{sample2,[a-z0-9.-]+}/done"
    params: bytes = lambda wildcards: metautils.getECS('foo','bytes','IsoModule'),
            mb = lambda wildcards: metautils.getECS('foo','mb','IsoModule'),
            gtf = "gencode.v28.annotation.gtf",
            reftx = "GRCh38_star",
            jobname = lambda wildcards: re.sub('\.','',wildcards.sample1+'_'+wildcards.sample2+'_sashimi')
    shell: """
            export untreated="{input.untreated}"
            export treated="{input.treated}"
            export manifest="{input.manifest}"
            export comparison="{wildcards.sample1}_vs_{wildcards.sample2}"
            export project="SRP091981-turbo"
            export bytes="{params.bytes}"
            export gtf="{params.gtf}" 
            export reftx="{params.reftx}"
            export gtf="{params.gtf}"
            export sample1="{wildcards.sample1}"
            export sample2="{wildcards.sample2}"
            export project="SRP091981-turbo"
            export destination="SRP091981-sashimi"
            scripts/sashimimanifest.sh
            """

rule bamindex:
    input: RAWDIR+"/{sample}.Aligned.sortedByCoord.out.bam",
    output: RAWDIR+"/{sample}.Aligned.sortedByCoord.out.bam.bai"
    params: bytes = lambda wildcards: metautils.getECS(wildcards.sample,'bytes','samtoolsindex'),
            mb = lambda wildcards: metautils.getECS(wildcards.sample,'mb','samtoolsindex')
    shell: """
            export sample="{wildcards.sample}"
            export project="SRP091981"
            export bytes="{params.bytes}"
            sh scripts/samtoolsindex.sh
            """

rule subsample:
    input: RAWDIR+"/{sample}.Aligned.sortedByCoord.out.bam",
    output: RAWDIR+"/{sample}.Aligned.sortedByCoord.out.subsampled.bam",RAWDIR+"/{sample}.Aligned.sortedByCoord.out.subsampled.bam.bai"
    params: bytes = lambda wildcards: metautils.getECS(wildcards.sample,'bytes','samtoolssubsample'),
            mb = lambda wildcards: metautils.getECS(wildcards.sample,'mb','samtoolssubsample'),
            subfraction = 0.001
    shell: """
            export sample="{wildcards.sample}"
            export project="SRP091981"
            export bytes="{params.bytes}"
            export subfraction="{params.subfraction}"
            sh scripts/subsample.sh
            """

rule rmatsiso:
    input: untreated=expand(PROJECT_BUCKET+"/"+RAWDIR+"/{sampleids}.Aligned.sortedByCoord.out.{ext}", ext=['bam','bam.bai'], sampleids=metautils.getRunsFromSampleName("Untreated HCT116")),
           treated=expand(PROJECT_BUCKET+"/"+RAWDIR+"/{sampleids}.Aligned.sortedByCoord.out.{ext}", ext=['bam','bam.bai'], sampleids=metautils.getRunsFromSampleName("0.5 uM T3 treated HCT116")),
           manifest="untreatedvslowdose.manifest.txt"
            
    output: expand(PROCESSDIR+"/{sampleids}.Aligned.sortedByCoord.out.bam.IsoExon", sampleids=metautils.getRunsFromSampleName("Untreated HCT116")),
            expand(PROCESSDIR+"/{sampleids}.Aligned.sortedByCoord.out.bam.IsoExon", sampleids=metautils.getRunsFromSampleName("0.5 uM T3 treated HCT116")),
            
    shell:
        """
        rMATS-ISO.py --in-gtf GRCh38_star/genes.gtf --in-bam {input.manifest} -o {PROCESSDIR}
        """
        

#rule starindex:
#    output: ""
#STAR --runMode genomeGenerate --runThreadN 8 --genomeDir ./ --genomeFastaFiles genome.fa --sjdbGTFfile genes.gtf --sjdbOverhang 100
# ensemblfetch homo_sapiens
# The command above retrieves the human genome sequence to a file called. Homo_sapiens.GRCh37.63.dna.toplevel.fa. You can calculate the GSNAP indexes for this file with command:

# gmap_build -d human -D $WRKDIR/gsnap_indexes Homo_sapiens.GRCh37.63.dna.toplevel.fa
# The option -d defines the name of the GSNAP indexes and the option -D defines the location, where the indexes will be stored.

# Single-end and pair-end alignment

# Once the indexing is ready you can carry out the alignment for singe-end reads with command like:

# gsnap -t 4 -d human -D $WRKDIR/gsnap_indexes query.fastq 
# In the case of paired-end reads you should have read pairs in two matching fastq files:

# gsnap -t 4 -d human -D $WRKDIR/gsnap_indexes query1.fastq query2.fastq 
# By default GSNAP uses its' own output format. To produce the alignment in SAM format please use option -A sam

# gsnap -t 4 -d human -D $WRKDIR/gsnap_indexes query1.fastq query2.fastq -A sam

#gmap_build -d g1k_v37 -k 15 -s none human_g1k_v37.fasta

# GSNAP_INDEX_DIR = "gsnap_index"
# GSNAP_INDEX_NAME = "hg37"
# rule gsnap_align:
#     input: pair1="{file}_1.fq",pair2="{file}_2.fq", gsnap = GSNAP_INDEX_DIR + "/" + GSNAP_INDEX_NAME
#     output: "{file}.bam"
#     threads: 4
#     shell:
#         """
#         gsnap -d g1k_v37 --gunzip -t {threads} -A sam -B 2 -N 1 {input.pair1} {input.pair2} | samtools view -bS - | samtools sort -  > {output} && samtools index {output}
#         """
        
#fix mates
# java -jar picard.jar FixMateInformation \
#       I=input.bam \
#       O=fixed_mate.bam \
#       ADD_MATE_CIGAR=true
       
       
rule map_mRNA:
    input: index="bacs_index", mrna="mRNA.fasta"
    output: "mRNA.sorted.bam"
    params: sge_opts="-l mfree=4G -pe serial 4"
    shell: "gmap -D `pwd` -d {input.index} -B 4 --min-identity=0.99 -t 4 --nofails --npaths=1 -A -f samse {input.mrna} 2> /dev/null | samtools view -Sbu -q 30 - | samtools sort -o {output} -O bam -T tmp_sort -"

rule build_GMAP_index:
    input: "clones.fasta"
    output: "bacs_index"
    params: sge_opts="-l mfree=5G"
    log: "gmap_build.log"
    shell: "gmap_build -D `pwd` -d {output} -k 15 -b 12 {input} &> {log}"
    
rule minimap_idx:
    input:
        config["genome"]["fasta"]
    output:
        config["genome"]["minimap_idx"]
    threads:
        config["minimap_idx"]["threads"]
    log:
        "logs/minimap_idx.log"
    benchmark:
        "benchmark/minimap_idx.benchmark.txt"
    params:
        minimap=config["exe_files"]["minimap2"]
    shell:
        "{params.minimap} -x splice {input} -d {output} -t {threads} 2> {log}"

rule sam_novel_gtf:
    input:
        sam=RAWDIR+"/{sample}.sam",
    output:
        filtered_bam=RAWDIR+"/{sample}.filtered.bam",
        sam_gtf=RAWDIR+"/{sample}_sam_novel.gtf"
    threads:
        config["novel_gtf"]["threads"]
    log:
        RAWDIR+"/{sample}.lr2rmats.log"
    # benchmark:
    #     "benchmark/{sample}.novel_gtf.benchmark.txt"
    params:
        aln_cov=config["lr2rmats"]["aln_cov"],
        iden_frac=config["lr2rmats"]["iden_frac"],
        sec_rat=config["lr2rmats"]["sec_rat"],
        mb = 10000
    shell:
        """
            export sample="{wildcards.sample}"
            export project=SRP091981" "aln_cov={params.aln_cov}" "iden_frac={params.iden_frac}" "sec_rat={params.sec_rat}" \
            --ebs /mnt/my-ebs:500:st1:ext4 \
            scripts/lr2rmats.sh && touch {output}
        """

