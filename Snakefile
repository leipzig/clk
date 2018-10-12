from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
from snakemake.utils import R
from utils import metautils
import yaml
import boto3
import re


#snakemake --no-shared-fs --default-remote-provider S3 --default-remote-prefix panorama-clk-repro  --cluster ./slurm_scheduler.py --cluster-status ./eric_status.py  -j 20 --cluster-config slurm_cluster_spec.yaml onemini

configfile: "config.yaml"

PROJECT_BUCKET = 'panorama-clk-repro'

def s3client():
    return boto3.client('s3',
        aws_access_key_id = os.environ['AWS_ACCESS_KEY_ID'],
        aws_secret_access_key = os.environ['AWS_SECRET_ACCESS_KEY']
    )
s3_boto2 = S3RemoteProvider()
s3_boto2.keep_local=True #Keep local copies of remote input files
#s3_boto2.stay_on_remote=True #don't download the file at all

LOCAL_SCRATCH = "/scratch"
RAWDIR="SRP091981"
PROCESSDIR="process"
SRAFILES = [line.rstrip() for line in open("SraAccList.txt")]
ILLUMINA_SRA = metautils.illuminaRuns()
PACBIO_SRA = metautils.pacbioRuns()

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

rule getapair:
    output: RAWDIR+"/{accession}_1.fastq.gz", RAWDIR+"/{accession}_2.fastq.gz"
    run:
        print("fastq-dump --split-3 --gzip {0} -O {1}".format(wildcards.accession,RAWDIR))
        shell("fastq-dump --split-3 --gzip {0} -O {1}".format(wildcards.accession,RAWDIR))
        #shell("fasterq-dump --split-3 {0} -O {1}".format(wildcards.accession,RAWDIR))
        #shell("gzip -1 {0}/{1}_1.fastq {0}/{1}_2.fastq".format(RAWDIR,wildcards.accession))

rule getasingleton:
    output: RAWDIR+"/{accessiondf -h }.fastq.gz"
    run:
        print("fastq-dump --split-3 --gzip {0} -O {1}".format(wildcards.accession,RAWDIR))
        shell("fastq-dump --split-3 --gzip {0} -O {1}".format(wildcards.accession,RAWDIR))

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

rule star_align:
    input: RAWDIR+"/{sample}_1.fastq.gz", RAWDIR+"/{sample}_2.fastq.gz"
    output: RAWDIR+"/{sample}.Aligned.sortedByCoord.out.bam",
            RAWDIR+"/{sample}.Log.final.out",
            RAWDIR+"/{sample}.Log.out",
            RAWDIR+"/{sample}.Log.progress.out",
            RAWDIR+"/{sample}.SJ.out.tab"
    params: bytes = lambda wildcards: metautils.getECS(wildcards.sample,'bytes','STAR'),
            mb = lambda wildcards: metautils.getECS(wildcards.sample,'mb','STAR')
    shell: """
            ./batchit submit \
            --image 977618787841.dkr.ecr.us-east-1.amazonaws.com/star:latest  \
            --role ecsTaskRole  \
            --queue when_you_get_to_it \
            --jobname star_{wildcards.sample} \
            --cpus 16 \
            --mem {params.mb} \
            --envvars "sample={wildcards.sample}" "project=SRP091981" "bytes={params.bytes}" \
            --ebs /mnt/my-ebs:500:st1:ext4 \
            star.sh
            touch {output}
            """
#> {wildcards.sample}.runid 2>&1 

rule helloworld:
    output: "helloworld.txt" #PROJECT_BUCKET+"/hello.txt"
    shell: """
            ./batchit submit \
            --image 977618787841.dkr.ecr.us-east-1.amazonaws.com/star:latest  \
            --role ecsTaskRole  \
            --queue when_you_get_to_it \
            --jobname hello \
            --cpus 1 \
            --mem 1000 \
            --ebs /mnt/my-ebs:500:st1:ext4 \
            hello.sh
            """

#output: PROJECT_BUCKET+"/hello_s3.txt")

rule hellos3:
    output: "hello_s3.txt"
    shell:
        """
        echo "hello world" > {output}
        aws s3 cp {output} s3://panorama-clk-repro/hello_s3.txt
        # echo "hello world" > hello_s3_tmp.txt
        # aws s3 cp hello_s3_tmp.txt s3://panorama-clk-repro/hello_s3.txt
        # touch {output}
        # echo 1234 >&1
        """
#cp hello_s3_tmp.txt {output}
#echo 1234 >&1
rule hellosimple:
    output: "hellosimple.txt"
    shell: """
            touch hellosimple.txt
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
        ./batchit submit \
        --image 977618787841.dkr.ecr.us-east-1.amazonaws.com/minimap2:latest  \
        --role ecsTaskRole  \
        --queue when_you_get_to_it \
        --jobname minimap_{wildcards.sample} \
        --cpus 16 \
        --mem {params.mb} \
        --envvars "sample={wildcards.sample}" "project=SRP091981" \
        --ebs /mnt/my-ebs:500:st1:ext4 \
        minimap.sh
        touch {output}
        """


# rule fetch_alignment:
#     input: s3_boto2.remote(PROJECT_BUCKET+"/"+RAWDIR+"/{sample}.Aligned.sortedByCoord.out.bam",keep_local=True)
#     output: RAWDIR+"/{sample}.Aligned.sortedByCoord.out.bam"
#     shell:
#         """
#         cp {input} {output}
#         """

#snakemake --no-shared-fs --default-remote-provider S3 --default-remote-prefix panorama-clk-repro  \
#panorama-clk-repro/SRP091981/untreated_vs_0.05.manifest.txt panorama-clk-repro/SRP091981/untreated_vs_0.1.manifest.txt panorama-clk-repro/SRP091981/untreated_vs_0.5.manifest.txt panorama-clk-repro/SRP091981/untreated_vs_1.0.manifest.txt
#panorama-clk-repro/SRP091981/untreated_vs_5.0.manifest.txt
#panorama-clk-repro/SRP091981/untreated184_vs_0.5-184.manifest.txt
#panorama-clk-repro/SRP091981/untreated184_vs_1.0-184.manifest.txt
#panorama-clk-repro/SRP091981/untreated184_vs_5.0-184.manifest.txt

rule anymanifest:
    output: manifest=RAWDIR+"/{sample1,[a-z0-9.-]+}_vs_{sample2,[a-z0-9.-]+}.manifest.txt"
    run:
        metautils.twoSampleComparisonManifest(metautils.dosageTable[wildcards.sample1],metautils.dosageTable[wildcards.sample2],output.manifest)


rule isomodule:
    input: bam=RAWDIR+"/{sample}.Aligned{ext}.bam", bai=RAWDIR+"/{sample}.Aligned{ext}.bam.bai"
    output: RAWDIR+"/{sample}.Aligned{ext}.bam.IsoExon", RAWDIR+"/{sample}.Aligned{ext}.bam.IsoMatrix"
    params: bytes = lambda wildcards: metautils.getECS(wildcards.sample,'bytes','IsoModule'),
            mb = lambda wildcards: metautils.getECS(wildcards.sample,'mb','IsoModule'),
            gtf = "level_1_protein_coding_genes.gtf"
    shell: """
            ./batchit submit \
            --image 977618787841.dkr.ecr.us-east-1.amazonaws.com/rmats-iso:latest  \
            --role ecsTaskRole  \
            --queue when_you_get_to_it \
            --jobname isomodule_{wildcards.sample} \
            --cpus 16 \
            --mem {params.mb} \
            --envvars "sample={wildcards.sample}.Aligned{wildcards.ext}" "project=SRP091981" "bytes={params.bytes}" "gtf={params.gtf}" \
            --ebs /mnt/my-ebs:500:st1:ext4 \
            isomodule.sh
            touch {output}
            """

# myoutput/
# myoutput/ISO_plot
# myoutput/ISO_module
# myoutput/ISO_module/PC3E-1.sorted.bam.IsoExon
# myoutput/ISO_module/GS689-3.sorted.bam.IsoMatrix
# myoutput/ISO_module/PC3E-2.sorted.bam.IsoMatrix
# myoutput/ISO_module/GS689-2.sorted.bam.IsoMatrix
# myoutput/ISO_module/PC3E-3.sorted.bam.IsoMatrix
# myoutput/ISO_module/GS689-1.sorted.bam.IsoMatrix
# myoutput/ISO_module/PC3E-1.sorted.bam.IsoMatrix
# myoutput/ISO_classify
# myoutput/ISO_classify/ISO_module_coor.txt
# myoutput/ISO_classify/ISO_module_gene.txt
# myoutput/ISO_classify/ISO_module_type_summary.txt
# myoutput/ISO_classify/ISO_module_type.txt
# myoutput/EM_out
# myoutput/EM_out/EM.out
#untreated=expand(RAWDIR+"/{sampleids}.Aligned.sortedByCoord.out.{ext}", sampleids=metautils.getRunsFromSampleName(lambda wildcards:metautils.getfulldosagename(wildcards.sample1)), ext=['bam','bam.bai']),
#treated=expand(RAWDIR+"/{sampleids}.Aligned.sortedByCoord.out.{ext}", sampleids=metautils.getRunsFromSampleName(lambda wildcards:metautils.getfulldosagename(wildcards.sample2)), ext=['bam','bam.bai']),

#snakemake --no-shared-fs --default-remote-provider S3 --default-remote-prefix panorama-clk-repro  --cluster ./slurm_scheduler.py --cluster-status ./eric_status.py  -j 50 --cluster-config slurm_cluster_spec.yaml \
#panorama-clk-repro/SRP091981/untreated_vs_0.1/done \
#panorama-clk-repro/SRP091981/untreated_vs_0.05/done \
#panorama-clk-repro/SRP091981/untreated_vs_1.0/done \
#panorama-clk-repro/SRP091981/untreated_vs_0.5/done \
#panorama-clk-repro/SRP091981/untreated_vs_5.0/done \
#panorama-clk-repro/SRP091981/untreated184_vs_0.5-184/done \
#panorama-clk-repro/SRP091981/untreated184_vs_1.0-184/done \
#panorama-clk-repro/SRP091981/untreated184_vs_5.0-184/done

rule isomodulefrommanifest:
    input: untreated=lambda wildcards: metautils.getBamsFromSampleName(metautils.getfulldosagename(wildcards.sample1),include_s3=RAWDIR),
           treated=lambda wildcards: metautils.getBamsFromSampleName(metautils.getfulldosagename(wildcards.sample2),include_s3=RAWDIR),
           manifest=RAWDIR+"/{sample1}_vs_{sample2}.manifest.txt"
    output: manifest=RAWDIR+"/{sample1}_vs_{sample2}/done"
    params: bytes = lambda wildcards: metautils.getECS('foo','bytes','IsoModule'),
            mb = lambda wildcards: metautils.getECS('foo','mb','IsoModule'),
            gtf = "gencode.v28.annotation.gtf",
            jobname = lambda wildcards: re.sub('\.','',wildcards.sample1+'_'+wildcards.sample2)
    shell: """
            ./batchit submit \
            --image 977618787841.dkr.ecr.us-east-1.amazonaws.com/rmats-source:latest  \
            --role ecsTaskRole  \
            --queue when_you_get_to_it \
            --jobname {params.jobname} \
            --cpus 16 \
            --mem {params.mb} \
            --envvars "untreated={input.untreated}" "treated={input.treated}" "manifest={input.manifest}" "comparison={wildcards.sample1}_vs_{wildcards.sample2}" "project=SRP091981" "bytes={params.bytes}" "gtf={params.gtf}" \
            --ebs /mnt/my-ebs:500:st1:ext4 \
            isomanifest.sh
            touch {output}
            """

rule rmatsmodule:
    input: untreated=expand(RAWDIR+"/{sampleids}.Aligned.sortedByCoord.out.bam.IsoMatrix", sampleids=metautils.getRunsFromSampleName("Untreated HCT116")),
           treated=expand(RAWDIR+"/{sampleids}.Aligned.sortedByCoord.out.bam.IsoExon", sampleids=metautils.getRunsFromSampleName("0.5 uM T3 treated HCT116"))

           
rule rmatsprep:
    input: untreated=s3_boto2.remote(expand(PROJECT_BUCKET+"/"+RAWDIR+"/{sampleids}.Aligned.sortedByCoord.out.bam", sampleids=metautils.getRunsFromSampleName("Untreated HCT116")),keep_local=True),
           treated=s3_boto2.remote(expand(PROJECT_BUCKET+"/"+RAWDIR+"/{sampleids}.Aligned.sortedByCoord.out.bam", sampleids=metautils.getRunsFromSampleName("0.5 uM T3 treated HCT116")),keep_local=True),
           manifest="untreatedvslowdose.manifest.txt"
    output: "prepped"
    shell: "ls -alt {input}; touch prepped"

    # input: untreated=s3_boto2.remote(expand(PROJECT_BUCKET+"/"+RAWDIR+"/{sampleids}.Aligned.sortedByCoord.out.bam", sampleids=metautils.getRunsFromSampleName("Untreated HCT116")),keep_local=True),
    #       treated=s3_boto2.remote(expand(PROJECT_BUCKET+"/"+RAWDIR+"/{sampleids}.Aligned.sortedByCoord.out.bam", sampleids=metautils.getRunsFromSampleName("0.5 uM T3 treated HCT116")),keep_local=True),
    
rule bamindex:
    input: RAWDIR+"/{sample}.Aligned.sortedByCoord.out.bam",
    output: RAWDIR+"/{sample}.Aligned.sortedByCoord.out.bam.bai"
    params: bytes = lambda wildcards: metautils.getECS(wildcards.sample,'bytes','samtoolsindex'),
            mb = lambda wildcards: metautils.getECS(wildcards.sample,'mb','samtoolsindex')
    shell: """
            ./batchit submit \
            --image 977618787841.dkr.ecr.us-east-1.amazonaws.com/samtools:latest  \
            --role ecsTaskRole  \
            --queue when_you_get_to_it \
            --jobname samtoolsindex_{wildcards.sample} \
            --cpus 1 \
            --mem {params.mb} \
            --envvars "sample={wildcards.sample}" "project=SRP091981" "bytes={params.bytes}" \
            --ebs /mnt/my-ebs:500:st1:ext4 \
            samtoolsindex.sh
            touch {output}
            """

rule subsample:
    input: RAWDIR+"/{sample}.Aligned.sortedByCoord.out.bam",
    output: RAWDIR+"/{sample}.Aligned.sortedByCoord.out.subsampled.bam",RAWDIR+"/{sample}.Aligned.sortedByCoord.out.subsampled.bam.bai"
    params: bytes = lambda wildcards: metautils.getECS(wildcards.sample,'bytes','samtoolssubsample'),
            mb = lambda wildcards: metautils.getECS(wildcards.sample,'mb','samtoolssubsample'),
            subfraction = 0.001
    shell: """
            ./batchit submit \
            --image 977618787841.dkr.ecr.us-east-1.amazonaws.com/samtools:latest  \
            --role ecsTaskRole  \
            --queue when_you_get_to_it \
            --jobname samtoolssubsample_{wildcards.sample} \
            --cpus 1 \
            --mem {params.mb} \
            --envvars "sample={wildcards.sample}" "project=SRP091981" "bytes={params.bytes}" "subfraction={params.subfraction}" \
            --ebs /mnt/my-ebs:500:st1:ext4 \
            subsample.sh
            touch {output}
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
            ./batchit submit \
            --image 977618787841.dkr.ecr.us-east-1.amazonaws.com/rmats-iso:latest  \
            --role ecsTaskRole  \
            --queue when_you_get_to_it \
            --jobname lr2rmats_{wildcards.sample} \
            --cpus 4 \
            --mem {params.mb} \
            --envvars "sample={wildcards.sample}" "project=SRP091981" "aln_cov={params.aln_cov}" "iden_frac={params.iden_frac}" "sec_rat={params.sec_rat}" \
            --ebs /mnt/my-ebs:500:st1:ext4 \
            lr2rmats.sh && touch {output}
        """

