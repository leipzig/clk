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
    input: expand(RAWDIR+"/{sampleids}.filtered.bam", sampleids=PACBIO_SRA)

rule pacbio_lr2rmats:
    input: expand(RAWDIR+"/{sampleids}_sam_novel.gtf", sampleids=PACBIO_SRA)

rule allfiles:
    input: ill = expand(RAWDIR+"/{sampleids}_{pair}.fastq.gz", sampleids=ILLUMINA_SRA, pair=[1,2]),
           pac = expand(RAWDIR+"/{sampleids}.fastq", sampleids=PACBIO_SRA)

rule fetchpair_from_aws:
    output: pair1 = RAWDIR+"/{accession}_1.fastq.gz",
            pair2 = RAWDIR+"/{accession}_2.fastq.gz"
    run:
        s3pair1 = metautils.st.loc[metautils.st['Run'] == wildcards.accession]['Pair1Filename'].to_string(index=False).replace(' ','')
        s3pair2 = metautils.st.loc[metautils.st['Run'] == wildcards.accession]['Pair2Filename'].to_string(index=False).replace(' ','')
        print(s3pair1)
        if wildcards.accession in metautils.pacbioRuns():
            raise ValueError('Not sure this should produce paired')
            #pair1 correspond to the bas.h5
            shell("wget https://clk-splicing.s3.amazonaws.com/SRP091981/{0}/{1}".format(wildcards.accession,s3pair1))
            shell("~/.local/bin/bash5tools.py {0} --outFilePrefix {1}".format(s3pair1,wildcards.accession))
        	#shell("bedtools bamtofastq -i {0} -fq {1} -fq2 {2}".format(output.pair1,output.pair2))
        else:
            shell("wget -O {2} https://clk-splicing.s3.amazonaws.com/SRP091981/{0}/{1}".format(wildcards.accession,s3pair1,output.pair1))
            shell("wget -O {2} https://clk-splicing.s3.amazonaws.com/SRP091981/{0}/{1}".format(wildcards.accession,s3pair2,output.pair2))

rule fetchpacbio_from_aws:
    #output: pair1 = RAWDIR+"/{accession}.fastq.gz"
    run:
        s3pair1 = metautils.st.loc[metautils.st['Run'] == wildcards.accession]['Pair1Filename'].to_string(index=False).replace(' ','')
        s3pair2 = metautils.st.loc[metautils.st['Run'] == wildcards.accession]['Pair2Filename'].to_string(index=False).replace(' ','')
        ena = metautils.st.loc[metautils.st['Run'] == wildcards.accession]['ena_fastq_http_1'].to_string(index=False).replace(' ','')
        print(s3pair1)
        if wildcards.accession in metautils.pacbioRuns():
            #pair1 correspond to the bas.h5
            shell("wget -O {0} {1}".format(output.pair1,ena))
            
            # terrible quality
            #shell("wget https://clk-splicing.s3.amazonaws.com/SRP091981/{0}/{1}".format(wildcards.accession,s3pair1))
            #shell("~/.local/bin/bash5tools.py {0} --outFilePrefix {1}/{2} --outType fastq".format(s3pair1,RAWDIR,wildcards.accession))


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
# STAR --runThreadN 16 --runMode genomeGenerate --genomeDir GRCh38_star_2.7.9 --genomeFastaFiles refs/GRCh38/Sequence/WholeGenomeFasta/genome.fa --sjdbGTFfile refs/GRCh38/Annotation/Genes.gencode/gencode.v38.annotation.gtf --sjdbOverhang 99
#/clk/refs/GRCh38/Annotation/Genes.gencode/gencode.v38.annotation.gtf
rule star_align:
    input: RAWDIR+"/{sample}_1.fastq.gz", RAWDIR+"/{sample}_2.fastq.gz"
    output: RAWDIR+"/{sample}.Aligned.sortedByCoord.out.bam",
            RAWDIR+"/{sample}.Log.final.out",
            RAWDIR+"/{sample}.Log.out",
            RAWDIR+"/{sample}.Log.progress.out",
            RAWDIR+"/{sample}.SJ.out.tab"
    threads: 16
    params: bytes = lambda wildcards: metautils.getECS(wildcards.sample,'bytes','STAR'),
            mb = lambda wildcards: metautils.getECS(wildcards.sample,'mb','STAR'),
            genomeDir = "GRCh38_star_2.7.9",
            gtf = "refs/GRCh38/Annotation/Genes.gencode/gencode.v38.annotation.gtf"
    shell: """
            STAR --runMode alignReads \
                --chimOutType WithinBAM \
                 --outSAMtype BAM SortedByCoordinate \
                 --limitBAMsortRAM {params.bytes} \
                 --readFilesCommand zcat \
                 --outFilterType BySJout   --outFilterMultimapNmax 20 \
                 --outFilterMismatchNmax 999   --alignIntronMin 25  \
                 --alignIntronMax 1000000   --alignMatesGapMax 1000000 \
                 --alignSJoverhangMin 8   --alignSJDBoverhangMin 5 \
                 --sjdbGTFfile {params.gtf} \
                 --genomeDir {params.genomeDir} \
                 --runThreadN {threads} \
                 --outFileNamePrefix SRP091981/{wildcards.sample}.  \
                 --readFilesIn  SRP091981/{wildcards.sample}_1.fastq.gz SRP091981/{wildcards.sample}_2.fastq.gz
            
            samtools index SRP091981/{wildcards.sample}.Aligned.sortedByCoord.out.bam
            """

# minimap mapping for long reads
rule minimap_map:
    input:
        "SRP091981/{sample}.fastq.gz"
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
        minimap2 -ax splice -ub -t {threads} GRCh38_minimap/genome.fa.smmi {input} > {output} 2> {wildcards.sample}.minimap.log
        """

rule generate_two_way_manifest:
    output: manifest=RAWDIR+"/{sample1,[a-z0-9.-]+}_vs_{sample2,[a-z0-9.-]+}.manifest.txt"
    run:
        metautils.twoSampleComparisonManifest(wildcards.sample1,wildcards.sample2,output.manifest,path_prefix="SRP091981")



rule run_rmatsiso_from_bam:
    input: bam=RAWDIR+"/{sample}.Aligned{ext}.bam", bai=RAWDIR+"/{sample}.Aligned{ext}.bam.bai"
    output: RAWDIR+"/{sample}.Aligned{ext}.bam.IsoExon", RAWDIR+"/{sample}.Aligned{ext}.bam.IsoMatrix"
    params: bytes = lambda wildcards: metautils.getECS(wildcards.sample,'bytes','IsoModule'),
            mb = lambda wildcards: metautils.getECS(wildcards.sample,'mb','IsoModule'),
            gtf = "level_1_protein_coding_genes.gtf"
    shell: """
            sample="{wildcards.sample}.Aligned{wildcards.ext}"
            project="SRP091981"
            bytes="{params.bytes}"
            gtf="{params.gtf}"
            sh scripts/isomodule.sh
            """

#gtf = "gencode.v28.annotation.gtf",

rule getbams:
    input: lambda wildcards: metautils.getBamsFromSampleName(wildcards.sample,path_prefix=RAWDIR)
    output: "{sample}.gotbams"
    shell: "touch {sample}.gotbams"

#metautils.getfulldosagename(wildcards.sample1) no longer necessary
rule run_rmatsiso_from_manifest:
    input: untreated=lambda wildcards: metautils.getBamsFromSampleName(wildcards.sample1,path_prefix=RAWDIR),
           treated=lambda wildcards: metautils.getBamsFromSampleName(wildcards.sample2,path_prefix=RAWDIR),
           manifest=RAWDIR+"/{sample1}_vs_{sample2}.manifest.txt"
    output: "results/iso_{sample1}_vs_{sample2}/ISO_module"
    params: bytes = lambda wildcards: metautils.getECS('foo','bytes','IsoModule'),
            mb = lambda wildcards: metautils.getECS('foo','mb','IsoModule'),
            gtf = "genes.gtf",
            jobname = lambda wildcards: re.sub('\.','',wildcards.sample1+'_'+wildcards.sample2),
            outdir = "results/iso_{sample1}_vs_{sample2}/"
    shell: """
            rMATS-ISO.py  --in-gtf GRCh38_star/{params.gtf} --in-bam {input.manifest} -o {params.outdir}
            """

rule run_rmatsturbo_from_manifest:
    input: untreated=lambda wildcards: metautils.getBamsFromSampleName(wildcards.sample1,path_prefix=RAWDIR),
           treated=lambda wildcards: metautils.getBamsFromSampleName(wildcards.sample2,path_prefix=RAWDIR),
           manifest=RAWDIR+"/{sample1}_vs_{sample2}.manifest.txt"
    output: RAWDIR+"-turbo/{sample1}_vs_{sample2}/rmats.out.txt"
    params: bytes = lambda wildcards: metautils.getECS('foo','bytes','IsoModule'),
            mb = lambda wildcards: metautils.getECS('foo','mb','IsoModule'),
            gtf = "refs/GRCh38/Annotation/Genes.gencode/genes.gtf",
            reftx = "GRCh38_star",
            jobname = lambda wildcards: re.sub('\.','',wildcards.sample1+'_'+wildcards.sample2),
            outdir = RAWDIR+"-turbo",
            tmpdir = "/tmp/{sample1}_vs_{sample2}"
    threads: 16
    shell: """
            rmats.py --readLength 150 --variable-read-length --allow-clipping --novelSS --nthread {threads} -t paired --gtf {params.gtf} --b1 <( python scripts/manifest_to_csl.py {input.manifest} 1 . )  \
            --b2 <( python scripts/manifest_to_csl.py {input.manifest} 2 . ) \
            --od {params.outdir}/{wildcards.sample1}_vs_{wildcards.sample2}  \
            --tmp {params.tmpdir} > {output}
            """

rule all_turbo:
    input: expand(RAWDIR+"-turbo/{sample1}_vs_{sample2}/rmats.out.txt",sample1="untreated",sample2=['0.05','0.1','0.5','1.0','treated'])
    
rule all_sashimi:
    input: expand(RAWDIR+"-sashimi/{sample1}_vs_{sample2}/done",sample1="untreated",sample2=['0.05','0.1','0.5','1.0','treated'])

#snakemake -j  SRP091981-sashimi/untreated_vs_treated/done panorama-clk-repro/SRP091981-sashimi/untreated_vs_0.1/done  panorama-clk-repro/SRP091981-sashimi/untreated_vs_0.5/done panorama-clk-repro/SRP091981-sashimi/untreated_vs_1.0/done 
rule run_rmatssashimi_from_manifest:
    input: untreated=lambda wildcards: metautils.getBamsFromSampleName(wildcards.sample1,path_prefix=RAWDIR),
           treated=lambda wildcards: metautils.getBamsFromSampleName(wildcards.sample2,path_prefix=RAWDIR),
           manifest=RAWDIR+"/{sample1}_vs_{sample2}.manifest.txt",
           rmats=RAWDIR+"-turbo/{sample1}_vs_{sample2}/rmats.out.txt"
    output: manifest=RAWDIR+"-sashimi/{sample1,[a-z0-9.-]+}_vs_{sample2,[a-z0-9.-]+}/done"
    params: bytes = lambda wildcards: metautils.getECS('foo','bytes','IsoModule'),
            mb = lambda wildcards: metautils.getECS('foo','mb','IsoModule'),
            gtf = "gencode.v28.annotation.gtf",
            reftx = "GRCh38_star",
            jobname = lambda wildcards: re.sub('\.','',wildcards.sample1+'_'+wildcards.sample2+'_sashimi')
    shell: """
            export EDITOR=emacs
            mkdir -p SRP091981-sashimi/{wildcards.sample1}_vs_{wildcards.sample2}/SE SRP091981-sashimi/{wildcards.sample1}_vs_{wildcards.sample2}/A5SS
            mkdir -p SRP091981-sashimi/{wildcards.sample1}_vs_{wildcards.sample2}/A3SS SRP091981-sashimi/{wildcards.sample1}_vs_{wildcards.sample2}/MXE
            mkdir -p SRP091981-sashimi/{wildcards.sample1}_vs_{wildcards.sample2}/RI
            Rscript scripts/filterRmats.R SRP091981-turbo/{wildcards.sample1}_vs_{wildcards.sample2}/SE.MATS.JC.txt omit
            Rscript scripts/filterRmats.R SRP091981-turbo/{wildcards.sample1}_vs_{wildcards.sample2}/A5SS.MATS.JC.txt omit
            Rscript scripts/filterRmats.R SRP091981-turbo/{wildcards.sample1}_vs_{wildcards.sample2}/A3SS.MATS.JC.txt omit
            Rscript scripts/filterRmats.R SRP091981-turbo/{wildcards.sample1}_vs_{wildcards.sample2}/MXE.MATS.JC.txt omit
            Rscript scripts/filterRmats.R SRP091981-turbo/{wildcards.sample1}_vs_{wildcards.sample2}/RI.MATS.JC.txt omit
            
            #some of these might fail
            set +e
            
            truncate -s 0 SRP091981-sashimi/{wildcards.sample1}_vs_{wildcards.sample2}/rmats-sashimi.out.txt
            
            #wc -l SRP091981-turbo/{wildcards.sample1}_vs_{wildcards.sample2}/*filtered.txt || echo "no eligible files"
            test -f SRP091981-turbo/{wildcards.sample1}_vs_{wildcards.sample2}/SE.MATS.JC.filtered.txt && rmats2sashimiplot --b1 `python scripts/manifest_to_csl.py {input.manifest} 1 .` \
                              --b2 `python scripts/manifest_to_csl.py {input.manifest} 2 .` \
                              -t SE \
                              -e SRP091981-turbo/SE.MATS.JC.filtered.txt \
                              --l1 {wildcards.sample1} \
                              --l2 {wildcards.sample2} \
                              --exon_s 1 \
                              --intron_s 5 \
                              -o SRP091981-sashimi/{wildcards.sample1}_vs_{wildcards.sample2}/SE >> SRP091981-sashimi/{wildcards.sample1}_vs_{wildcards.sample2}/rmats-sashimi.out.txt 2>&1
            
            test -f SRP091981-turbo/{wildcards.sample1}_vs_{wildcards.sample2}/A5SS.MATS.JC.filtered.txt && rmats2sashimiplot --b1 `python scripts/manifest_to_csl.py {input.manifest} 1 .` \
                              --b2 `python scripts/manifest_to_csl.py {input.manifest} 2 .` \
                              -t A5SS \
                              -e SRP091981-turbo/{wildcards.sample1}_vs_{wildcards.sample2}/A5SS.MATS.JC.filtered.txt \
                              --l1 {wildcards.sample1} \
                              --l2 {wildcards.sample2} \
                              --exon_s 1 \
                              --intron_s 5 \
                              -o SRP091981-sashimi/{wildcards.sample1}_vs_{wildcards.sample2}/A5SS > SRP091981-sashimi/{wildcards.sample1}_vs_{wildcards.sample2}/rmats-sashimi.out.txt 2>&1
            
            test -f SRP091981-turbo/A3SS.MATS.JC.filtered.txt && rmats2sashimiplot --b1 `python scripts/manifest_to_csl.py {input.manifest} 1 .` \
                              --b2 `python scripts/manifest_to_csl.py {input.manifest} 2 .` \
                              -t A3SS \
                              -e SRP091981-turbo/{wildcards.sample1}_vs_{wildcards.sample2}/A3SS.MATS.JC.filtered.txt \
                              --l1 {wildcards.sample1} \
                              --l2 {wildcards.sample2} \
                              --exon_s 1 \
                              --intron_s 5 \
                              -o SRP091981-sashimi/{wildcards.sample1}_vs_{wildcards.sample2}/A3SS >> SRP091981-sashimi/{wildcards.sample1}_vs_{wildcards.sample2}/rmats-sashimi.out.txt 2>&1
            
            test -f SRP091981-turbo/MXE.MATS.JC.filtered.txt && rmats2sashimiplot --b1 `python scripts/manifest_to_csl.py {input.manifest} 1 .` \
                              --b2 `python scripts/manifest_to_csl.py {input.manifest} 2 .` \
                              -t MXE \
                              -e SRP091981-turbo/{wildcards.sample1}_vs_{wildcards.sample2}/MXE.MATS.JC.filtered.txt \
                              --l1 {wildcards.sample1} \
                              --l2 {wildcards.sample2} \
                              --exon_s 1 \
                              --intron_s 5 \
                              -o SRP091981-sashimi/{wildcards.sample1}_vs_{wildcards.sample2}/MXE >> SRP091981-sashimi/{wildcards.sample1}_vs_{wildcards.sample2}/rmats-sashimi.out.txt 2>&1
            
            test -f SRP091981-turbo/RI.MATS.JC.filtered.txt && rmats2sashimiplot --b1 `python scripts/manifest_to_csl.py {input.manifest} 1 .` \
                              --b2 `python scripts/manifest_to_csl.py {input.manifest} 2 .` \
                              -t RI \
                              -e SRP091981-turbo/{wildcards.sample1}_vs_{wildcards.sample2}/RI.MATS.JC.filtered.txt \
                              --l1 {wildcards.sample1} \
                              --l2 {wildcards.sample2} \
                              --exon_s 1 \
                              --intron_s 5 \
                              -o SRP091981-sashimi/{wildcards.sample1}_vs_{wildcards.sample2}/RI >> SRP091981-sashimi/{wildcards.sample1}_vs_{wildcards.sample2}/rmats-sashimi.out.txt 2>&1
                              
            touch {output}
            """

rule bamindex:
    input: RAWDIR+"/{sample}.Aligned.sortedByCoord.out.bam",
    output: RAWDIR+"/{sample}.Aligned.sortedByCoord.out.bam.bai"
    params: bytes = lambda wildcards: metautils.getECS(wildcards.sample,'bytes','samtoolsindex'),
            mb = lambda wildcards: metautils.getECS(wildcards.sample,'mb','samtoolsindex')
    shell: """
            samtools index {input}
            """

rule subsample:
    input: RAWDIR+"/{sample}.Aligned.sortedByCoord.out.bam",
    output: RAWDIR+"/{sample}.Aligned.sortedByCoord.out.subsampled.bam",RAWDIR+"/{sample}.Aligned.sortedByCoord.out.subsampled.bam.bai"
    params: bytes = lambda wildcards: metautils.getECS(wildcards.sample,'bytes','samtoolssubsample'),
            mb = lambda wildcards: metautils.getECS(wildcards.sample,'mb','samtoolssubsample'),
            subfraction = 0.001
    shell: """
            sample="{wildcards.sample}"
            project="SRP091981"
            bytes="{params.bytes}"
            subfraction="{params.subfraction}"
            sh scripts/subsample.sh
            """

rule rmatsisooneoff:
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

rule rrnagtf:
    output: "GRCh38_star/rRNA_tx.gtf"
    shell: "curl https://raw.githubusercontent.com/zxl124/rRNA_gtfs/master/UCSC/hg38.gtf  > {output}"
    
rule sam_novel_gtf:
    input:
        sam=RAWDIR+"/{sample}.sam",
        gtf = "GRCh38_star/rRNA_tx.gtf"
    output:
        filtered_bam=RAWDIR+"/{sample}.filtered.bam",filtered_bai=RAWDIR+"/{sample}.filtered.bam.bai",
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
        lr2rmats filter {input.sam} -r GRCh38_star/rRNA_tx.gtf -v {params.aln_cov} -q {params.iden_frac} -s {params.sec_rat}  | samtools sort -@ {threads} > {output.filtered_bam}
        samtools index {output.filtered_bam}
        lr2rmats update-gtf {output.filtered_bam} GRCh38_star/genes.gtf > {output.sam_gtf}
        """

rule fetchtx:
    output: "refs/GRCh38/Sequence/Transcriptome/gencode.v38.transcripts.fa.gz"
    shell: " wget -O {output} http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.transcripts.fa.gz"

rule fetchgencodegtf:
    shell: "wget -O refs/GRCh38/Annotation/Genes.gencode/gencode.v38.annotation.gtf.gz http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz"

#some of the read lengths in SRP091981 are very short (<45)
# a lower k-mer index was chosen
rule salmonindex:
    input: "refs/GRCh38/Sequence/Transcriptome/gencode.v38.transcripts.fa.gz"
    output: "refs/gencode.salmon.v38/versionInfo.json"
    params: refdir = "refs/gencode.salmon.v38/"
    shell:
        """
        salmon index --gencode -t {input} -i {params.refdir} --kmerLen 17
        """
        
rule salmonquant:
    input: ref = "refs/gencode.salmon.v38/versionInfo.json", pair1 = RAWDIR+"/{sample}_1.fastq.gz", pair2 = RAWDIR+"/{sample}_2.fastq.gz"
    output: RAWDIR+"/salmon/paired/{sample}/quant.sf"
    params: outdir = RAWDIR+"/salmon/paired/{sample}"
    threads: 4
    shell:
        """
        salmon quant -i refs/gencode.salmon.v38 -l A --gcBias -1 {input.pair1} -2 {input.pair2} -p {threads} -o {params.outdir}
        """

rule salmonquantsingle:
    input: ref = "refs/gencode.salmon.v38/versionInfo.json", single = RAWDIR+"/{sample}.fastq.gz"
    output: RAWDIR+"/salmon/single/{sample}/quant.sf"
    params: outdir = RAWDIR+"/salmon/single/{sample}"
    threads: 4
    shell:
        """
        salmon quant -i refs/gencode.salmon.v38 -l A --minAssignedFrags 1 -r {input.single} -p {threads} -o {params.outdir}
        """
        
rule salmon:
    input: expand(RAWDIR+"/salmon/paired/{sampleids}/quant.sf", sampleids=ILLUMINA_SRA), expand(RAWDIR+"/salmon/single/{sampleids}/quant.sf", sampleids=PACBIO_SRA) 
    

rule suppaindex:
    # It requires python3 and pandas library (pip install pandas)
    # -k indicates the row used as the index
    # -f indicates the column to be extracted from the Salmon output
    #cat iso_tpm_deflinehell.txt | sed -e 's/|\S*//' > iso_tpm.txt
    shell: """
           multipleFieldSelection.py -i SRP091981/salmon/*/*/quant.sf -k 1 -f 4 -o iso_tpm.txt
           """

rule suppapsi:
    input: "iso_tpm.txt"
    output: "results/psiPerIsoform_isoform.psi"
    shell: "suppa.py psiPerIsoform -g refs/GRCh38/Annotation/Genes.gencode/gencode.v38.annotation.gtf -e {input} -o results/psiPerIsoform"

rule arribarefs:
    shell:
        """
        /home/ec2-user/miniconda3/envs/clk/var/lib/arriba/download_references.sh GRCh38+GENCODE28
        """

rule callFusionsArriba:
    input: RAWDIR+"/{sample}.Aligned.sortedByCoord.out.bam"
    output: RAWDIR+"/{sample}.fusions.tsv"
    shell:
        """
        ARRIBA_FILES=$CONDA_PREFIX/var/lib/arriba
        arriba -x {input} \
           -g GRCh38_star/genes.gtf -a GRCh38.fa \
           -b $ARRIBA_FILES/blacklist_hg38_GRCh38_v2.1.0.tsv.gz -k $ARRIBA_FILES/known_fusions_hg38_GRCh38_v2.1.0.tsv.gz \
           -p $ARRIBA_FILES/protein_domains_hg38_GRCh38_v2.1.0.gff3 \
           -o {output}
        """

rule fusioned:
    input: expand(RAWDIR+"/{sampleids}.fusions.tsv", sampleids=ILLUMINA_SRA)