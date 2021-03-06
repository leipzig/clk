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

rule all:
    input: "results/fusionCount.txt", "SRP091981-sashimi/untreated_vs_treated/done"

rule onemini:
    input: RAWDIR+"/SRR5009429_sam_novel.gtf", RAWDIR+"/SRR5009429.lr2rmats.log", RAWDIR+"/SRR5009429.filtered.bam"

rule justone:
    input: RAWDIR+"/SRR5009515.Aligned.sortedByCoord.out.md.bam"

rule bigboy:
    input: RAWDIR+"/SRR5009517.Aligned.sortedByCoord.out.md.bam"

   
rule s3_illumina_files:
    input: expand(RAWDIR+"/{sampleids}_{pair}.fastq.gz", sampleids=ILLUMINA_SRA, pair=[1,2])

rule illumina_align:
    input: expand(RAWDIR+"/{sampleids}.Aligned.sortedByCoord.out.md.bam", sampleids=ILLUMINA_SRA)

rule illumina_index:
    input: expand(RAWDIR+"/{sampleids}.Aligned.sortedByCoord.out.md.bam.bai", sampleids=ILLUMINA_SRA)


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
    output: RAWDIR+"/{sample}.Aligned.sortedByCoord.out.md.bam",
            RAWDIR+"/{sample}.Log.final.out",
            RAWDIR+"/{sample}.Log.out",
            RAWDIR+"/{sample}.Log.progress.out",
            RAWDIR+"/{sample}.SJ.out.tab"
    threads: 1
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
                 --outFileNamePrefix SRP091981/{wildcards.sample}.nochim.  \
                 --readFilesIn  SRP091981/{wildcards.sample}_1.fastq.gz SRP091981/{wildcards.sample}_2.fastq.gz
            
            samtools index SRP091981/{wildcards.sample}.Aligned.sortedByCoord.out.md.bam
            """
            
rule chimera_star_align:
    input: RAWDIR+"/{sample}_1.fastq.gz", RAWDIR+"/{sample}_2.fastq.gz"
    output: RAWDIR+"/{sample}.chim.Aligned.sortedByCoord.out.md.bam",
            RAWDIR+"/{sample}.chim.Log.final.out",
            RAWDIR+"/{sample}.chim.Log.out",
            RAWDIR+"/{sample}.chim.Log.progress.out",
            RAWDIR+"/{sample}.chim.SJ.out.tab"
    threads: 16
    params: bytes = lambda wildcards: metautils.getECS(wildcards.sample,'bytes','STAR'),
            mb = lambda wildcards: metautils.getECS(wildcards.sample,'mb','STAR'),
            genomeDir = "GRCh38_star_2.7.9",
            gtf = "refs/GRCh38/Annotation/Genes.gencode/gencode.v38.annotation.gtf"
    shell: """
            STAR --runMode alignReads \
            --outFilterMultimapNmax 50 \
            --peOverlapNbasesMin 10 \
            --alignSplicedMateMapLminOverLmate 0.5 \
            --alignSJstitchMismatchNmax 5 -1 5 5 \
            --chimSegmentMin 10 \
            --chimOutType WithinBAM HardClip \
            --chimJunctionOverhangMin 10 \
            --chimScoreDropMax 30 \
            --chimScoreJunctionNonGTAG 0 \
            --chimScoreSeparation 1 \
            --chimSegmentReadGapMax 3 \
            --chimMultimapNmax 50 \
            --genomeDir {params.genomeDir} \
            --runThreadN {threads} \
            --outFileNamePrefix SRP091981/{wildcards.sample}.chim.  \
            --readFilesIn  SRP091981/{wildcards.sample}_1.fastq.gz SRP091981/{wildcards.sample}_2.fastq.gz
            
            samtools index SRP091981/{wildcards.sample}.chim.Aligned.sortedByCoord.out.md.bam
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

rule iso_classify:
    input: isoexon="results/iso_untreated_vs_{dose}/ISO_module/SRR5009496.Aligned.sortedByCoord.out.md.bam.IsoExon",
            emout="results/iso_untreated_vs_{dose}/EM_out/EM.out"
    output: summary="results/iso_untreated_vs_{dose}/ISO_classify/ISO_module_type_summary.txt",
            stype="results/iso_untreated_vs_{dose}/ISO_classify/ISO_module_type.txt",
            coor="results/iso_untreated_vs_{dose}/ISO_classify/ISO_module_coor.txt",
            gene="results/iso_untreated_vs_{dose}/ISO_classify/ISO_module_gene.txt",
    shell:
       """
       mkdir -p results/iso_untreated_vs_0.5/ISO_classify/
       python2.7 rMATS-ISO-master/ISOClassify/IsoClass.py {input.isoexon} {output.summary} {output.stype}
       python2.7 rMATS-ISO-master/ISOPlot/IsoPlot.py {input.emout} {input.isoexon} {output.coor} {output.gene}
       """

#gtf = "gencode.v28.annotation.gtf",

rule getbams:
    input: lambda wildcards: metautils.getBamsFromSampleName(wildcards.sample,path_prefix=RAWDIR)
    output: "{sample}.gotbams"
    shell: "touch {sample}.gotbams"

#metautils.getfulldosagename(wildcards.sample1) no longer necessary


#results/iso_untreated_vs_0.05/ISO_module/done results/iso_untreated_vs_0.5/ISO_module/done results/iso_untreated_vs_5.0/ISO_module/done results/iso_untreated_vs_0.1/ISO_module/done results/iso_untreated_vs_1.0/ISO_module/done results/iso_untreated_vs_treated/ISO_module/done

rule run_rmatsiso_from_manifest:
    input: untreated=lambda wildcards: metautils.getBamsFromSampleName(wildcards.sample1,path_prefix=RAWDIR,platform='ILLUMINA'),
           treated=lambda wildcards: metautils.getBamsFromSampleName(wildcards.sample2,path_prefix=RAWDIR,platform='ILLUMINA'),
           manifest=RAWDIR+"/{sample1}_vs_{sample2}.manifest.txt"
    output: "results/iso_{sample1}_vs_{sample2}/ISO_module/done"
    params: bytes = lambda wildcards: metautils.getECS('foo','bytes','IsoModule'),
            mb = lambda wildcards: metautils.getECS('foo','mb','IsoModule'),
            gtf = "refs/Homo_sapiens.GRCh38.104.chred.gtf",
            jobname = lambda wildcards: re.sub('\.','',wildcards.sample1+'_'+wildcards.sample2),
            outdir = "results/iso_{sample1}_vs_{sample2}/"
    shell: """
            python rMATS-ISO-master/rMATS-ISO.py  module  --gtf {params.gtf} --bam {input.manifest} -o {params.outdir}
            touch {output}
            """

rule run_rmatsiso_stat:
    input: module="results/iso_{sample1}_vs_{sample2}/ISO_module",
           manifest=RAWDIR+"/{sample1}_vs_{sample2}.manifest.txt"
    output: "results/iso_{sample1}_vs_{sample2}/EM_out/EM.out"
    params: outdir = "results/iso_{sample1}_vs_{sample2}/"
    threads: 4
    shell:
            """
            python rMATS-ISO-master/rMATS-ISO.py stat --bam {input.manifest} -o {params.outdir}
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
    input: RAWDIR+"/{sample}.Aligned.sortedByCoord.out.md.bam",
    output: RAWDIR+"/{sample}.Aligned.sortedByCoord.out.md.bam.bai"
    params: bytes = lambda wildcards: metautils.getECS(wildcards.sample,'bytes','samtoolsindex'),
            mb = lambda wildcards: metautils.getECS(wildcards.sample,'mb','samtoolsindex')
    shell: """
            samtools index {input}
            """

rule subsample:
    input: RAWDIR+"/{sample}.Aligned.sortedByCoord.out.md.bam",
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
    input: untreated=expand(RAWDIR+"/{sampleids}.Aligned.sortedByCoord.out.md.{ext}", ext=['bam','bam.bai'], sampleids=metautils.getRunsFromSampleName("Untreated HCT116")),
           treated=expand(RAWDIR+"/{sampleids}.Aligned.sortedByCoord.out.md.{ext}", ext=['bam','bam.bai'], sampleids=metautils.getRunsFromSampleName("0.5 uM T3 treated HCT116")),
           manifest="SRP091981/untreatedvslowdose.manifest.txt"
            
    output: expand(PROCESSDIR+"/{sampleids}.Aligned.sortedByCoord.out.md.bam.IsoExon", sampleids=metautils.getRunsFromSampleName("Untreated HCT116")),
            expand(PROCESSDIR+"/{sampleids}.Aligned.sortedByCoord.out.md.bam.IsoExon", sampleids=metautils.getRunsFromSampleName("0.5 uM T3 treated HCT116")),
            
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

rule suppaioe:
    output: expand("clk_{event}_strict.ioe",event=['A5','AF','AL','A3','MX','RI','SE'])
    shell:
        """
        suppa.py generateEvents -i   refs/GRCh38/Annotation/Genes.gencode/gencode.v38.annotation.gtf -o clk -f ioe -e SE SS MX RI FL
        """

#suppa.py psiPerEvent --ioe-file clk.ioe_A3_strict.ioe --expression-file iso_tpm.txt -o clk.A3_strict.psiperevent

#suppa.py psiPerEvent --ioe-file clk.ioe_A3_strict.ioe --expression-file iso_tpm.txt -o clk.A3_strict

rule makepsiperevent:
    input: tpm="iso_tpm.txt", ioe="clk_{event}_strict.ioe"
    output: "clk_{event}_strict.psi"
    params: core="clk_{event}_strict"
    shell:
        """
        suppa.py psiPerEvent --ioe-file {input.ioe} --expression-file {input.tpm} -o {params.core}
        """



rule isosplitall:
    input: "{file}.{ext}"
    output: expand("{seqtype}_{{file}}_{dose}.{{ext}}",seqtype=['all','ill'],dose=['untreated','05','10','5','005','1','01'])
    shell: """
./scripts/subsetcols.py {wildcards.file}.{wildcards.ext} 'SRR5009494,SRR5009490,SRR5009482,SRR5009470,SRR5009451,SRR5009505,SRR5009501,SRR5009496,SRR5009387,SRR5009390,SRR5009396,SRR5009398,SRR5009434,SRR5009429,SRR5009425,SRR5009404,SRR5009463,SRR5009465,SRR5009468,SRR5009438,SRR5009534,SRR5009521,SRR5009474'  all_{wildcards.file}_untreated.{wildcards.ext}
./scripts/subsetcols.py {wildcards.file}.{wildcards.ext} 'SRR5009526,SRR5009381'  all_{wildcards.file}_005.{wildcards.ext}
./scripts/subsetcols.py {wildcards.file}.{wildcards.ext} 'SRR5009378,SRR5009464'  all_{wildcards.file}_01.{wildcards.ext}
./scripts/subsetcols.py {wildcards.file}.{wildcards.ext} 'SRR5009513,SRR5009519,SRR5009528,SRR5009371,SRR5009373,SRR5009375,SRR5009376,SRR5009379,SRR5009405,SRR5009409,SRR5009412,SRR5009403,SRR5009423,SRR5009424,SRR5009435,SRR5009442,SRR5009448,SRR5009469,SRR5009487,SRR5009492,SRR5009515,SRR5009516,SRR5009377,SRR5009459'  all_{wildcards.file}_05.{wildcards.ext}
./scripts/subsetcols.py {wildcards.file}.{wildcards.ext} 'SRR5009514,SRR5009491,SRR5009462,SRR5009453'  all_{wildcards.file}_1.{wildcards.ext}
./scripts/subsetcols.py {wildcards.file}.{wildcards.ext} 'SRR5009504,SRR5009440,SRR5009481,SRR5009447,SRR5009436,SRR5009460,SRR5009428,SRR5009432,SRR5009433,SRR5009414,SRR5009416,SRR5009394,SRR5009399,SRR5009402,SRR5009380,SRR5009384,SRR5009523,SRR5009530,SRR5009532,SRR5009437,SRR5009392,SRR5009444,SRR5009479'  all_{wildcards.file}_5.{wildcards.ext}
./scripts/subsetcols.py {wildcards.file}.{wildcards.ext} 'SRR5009509,SRR5009383'  all_{wildcards.file}_10.{wildcards.ext}

./scripts/subsetcols.py {wildcards.file}.{wildcards.ext} 'SRR5009496,SRR5009521,SRR5009474'  ill_{wildcards.file}_untreated.{wildcards.ext}
./scripts/subsetcols.py {wildcards.file}.{wildcards.ext} 'SRR5009526,SRR5009381'  ill_{wildcards.file}_005.{wildcards.ext}
./scripts/subsetcols.py {wildcards.file}.{wildcards.ext} 'SRR5009378,SRR5009464'  ill_{wildcards.file}_01.{wildcards.ext}
./scripts/subsetcols.py {wildcards.file}.{wildcards.ext} 'SRR5009487,SRR5009515,SRR5009377,SRR5009459'  ill_{wildcards.file}_05.{wildcards.ext}
./scripts/subsetcols.py {wildcards.file}.{wildcards.ext} 'SRR5009514,SRR5009491,SRR5009462,SRR5009453'  ill_{wildcards.file}_1.{wildcards.ext}
./scripts/subsetcols.py {wildcards.file}.{wildcards.ext} 'SRR5009437,SRR5009392,SRR5009444,SRR5009479'  ill_{wildcards.file}_5.{wildcards.ext}
./scripts/subsetcols.py {wildcards.file}.{wildcards.ext} 'SRR5009509,SRR5009383'  ill_{wildcards.file}_10.{wildcards.ext} 
            """


rule isotargets:
    input: "iso_tpm_10.txt","clk.A3_strict_10.psiperevent"

rule diffsplice:
    input: ioe="clk.ioe_{sptype}_strict.ioe", c1psi="ill_clk_{sptype}_strict_{dose}.psi", utpsi="ill_clk_{sptype}_strict_untreated.psi", c1tpm="ill_iso_tpm_{dose}.txt", uttpm="ill_iso_tpm_untreated.txt"
    output: "diff_{sptype}_strict_{dose}.dpsi","diff_{sptype}_strict_{dose}.psivec"
    params: core = "diff_{sptype}_strict_{dose}"
    shell: """
            suppa.py diffSplice -m classical --input {input.ioe} --psi {input.utpsi} {input.c1psi}  --tpm {input.uttpm} {input.c1tpm}  --area 1000 --lower-bound 0.0 -gc -o {params.core}
            """

rule alldiff:
    input: expand("diff_{event}_strict_{dose}.dpsi",event=['A5','AF','AL','A3','MX','RI','SE'],dose=['05','10','5','005','1','01'])

rule allsimplified:
    input: expand("diff_{event}_strict_{dose}.dpsi.simplified",event=['A5','AF','AL','A3','MX','RI','SE'],dose=['05','10','5','005','1','01'])
    output: "all.simplified.csv"
    shell: "csvstack -t diff*simplified > all.simplified.csv"

rule simplifydiff:
    input: "diff_{event}_strict_{dose}.dpsi"
    output: "diff_{event}_strict_{dose}.dpsi.simplified"
    shell: "./scripts/simplifydiff.py {input} {wildcards.event} {wildcards.dose} {output}"

rule suppapsi:
    input: "iso_tpm.txt"
    output: "results/psiPerIsoform_isoform.psi"
    shell: "suppa.py psiPerIsoform -g refs/GRCh38/Annotation/Genes.gencode/gencode.v38.annotation.gtf -e {input} -o results/psiPerIsoform"

#hg38 is chr, goes with gencode used to regenerate STAR indexes. sigh.
rule arribarefs:
    shell:
        """
        /home/ec2-user/miniconda3/envs/clk/var/lib/arriba/download_references.sh GRCh38+GENCODE28
        /home/ec2-user/miniconda3/envs/clk/var/lib/arriba/download_references.sh hg38+GENCODE28
        """

#es refs/GRCh38/Sequence/WholeGenomeFasta/genome.fa --sjdbGTFfile refs/GRCh38/Annotation/Genes.gencode/gencode.v38.annotation.gtf --sjdbOverhang 99
#/clk/refs/GRCh38/Annotation/Genes.gencode/gencode.v38.annotation.gtf
#/home/ec2-user/miniconda3/envs/clk/var/lib/arriba/blacklist_hg38_h38_v2.1.0.tsv.gz chrified
#/home/ec2-user/miniconda3/envs/clk/var/lib/arriba/protein_domains_hg38_h38_v2.1.0.gff3 
rule callFusionsArriba:
    input: pair1 = RAWDIR+"/{sample}_1.fastq.gz", pair2 =  RAWDIR+"/{sample}_2.fastq.gz"
    output: RAWDIR+"/{sample}.fusions.tsv"
    threads: 8
    shell:
        """
        ARRIBA_FILES=$CONDA_PREFIX/var/lib/arriba
        STAR \
        --runThreadN 8 \
        --outTmpDir {wildcards.sample}.tmp \
        --genomeDir STAR_index_hg38_GENCODE28 --genomeLoad NoSharedMemory \
        --readFilesIn {input}  --readFilesCommand zcat \
        --outStd BAM_Unsorted --outSAMtype BAM Unsorted --outSAMunmapped Within --outBAMcompression 0 \
        --outFilterMultimapNmax 50 --peOverlapNbasesMin 10 --alignSplicedMateMapLminOverLmate 0.5 --alignSJstitchMismatchNmax 5 -1 5 5 \
        --chimSegmentMin 10 --chimOutType WithinBAM HardClip --chimJunctionOverhangMin 10 --chimScoreDropMax 30 \
        --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 --chimSegmentReadGapMax 3 --chimMultimapNmax 50 | \
        arriba -x /dev/stdin \
           -g GENCODE28.gtf  -a hg38.fa \
           -b $ARRIBA_FILES/blacklist_hg38_h38_v2.1.0.tsv.gz -k $ARRIBA_FILES/known_fusions_hg38_h38_v2.1.0.tsv.gz \
           -p $ARRIBA_FILES/protein_domains_hg38_h38_v2.1.0.gff3 \
           -o {output}
        """

rule callFusionsArribaPacbio:
    input: pair1 = RAWDIR+"/{sample}.fastq.gz"
    output: RAWDIR+"/{sample,[^_]+}.fusions.pb.tsv"
    threads: 8
    shell:
        """
        ARRIBA_FILES=$CONDA_PREFIX/var/lib/arriba
        STAR \
        --runThreadN 8 \
        --outTmpDir {wildcards.sample}.tmp \
        --genomeDir STAR_index_hg38_GENCODE28 --genomeLoad NoSharedMemory \
        --readFilesIn {input} --readFilesCommand zcat \
        --outStd BAM_Unsorted --outSAMtype BAM Unsorted --outSAMunmapped Within --outBAMcompression 0 \
        --outFilterMultimapNmax 50 --peOverlapNbasesMin 10 --alignSplicedMateMapLminOverLmate 0.5 --alignSJstitchMismatchNmax 5 -1 5 5 \
        --chimSegmentMin 10 --chimOutType WithinBAM HardClip --chimJunctionOverhangMin 10 --chimScoreDropMax 30 \
        --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 --chimSegmentReadGapMax 3 --chimMultimapNmax 50 | \
        arriba -x /dev/stdin \
           -g GENCODE28.gtf  -a hg38.fa \
           -b $ARRIBA_FILES/blacklist_hg38_h38_v2.1.0.tsv.gz -k $ARRIBA_FILES/known_fusions_hg38_h38_v2.1.0.tsv.gz \
           -p $ARRIBA_FILES/protein_domains_hg38_h38_v2.1.0.gff3 \
           -o {output}
        """
        
rule fusioned:
    input: expand(RAWDIR+"/{sampleids}.fusions.tsv", sampleids=ILLUMINA_SRA)
    output: "results/fusionCounts.txt"
    shell:
        """
        echo "Type\tRun\tcount" > results/fusionCounts.txt
        grep -c '^[^#]' SRP091981/*.fusions.tsv | sed -e 's/SRP091981\//all\t/' | sed -e 's/.fusions.tsv:/\t/' >> results/fusionCounts.txt
        grep -c 'read-through' SRP091981/*.fusions.tsv | sed -e 's/SRP091981\//readthrough\t/' | sed -e 's/.fusions.tsv:/\t/' >> results/fusionCounts.txt
        """
