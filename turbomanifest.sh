set -euo pipefail
cd $TMPDIR

echo "downloading bams.."
for f in ${untreated}; do aws s3 cp s3://$f .; done;
for f in ${treated}; do aws s3 cp s3://$f .; done;

aws s3 cp s3://${manifest} manifest.list

echo "I am here:"
pwd
ls -alt .

echo "downloading gtf.."
mkdir GRCh38_star
aws s3 cp s3://panorama-refs/GRCh38_star/${gtf} GRCh38_star/

    # usage: python rmats.py [options] arg1 arg2
	
    # optional arguments:
		  #-h, --help            show this help message and exit
		  #--version             Version.
		  #--gtf GTF             An annotation of genes and transcripts in GTF format.
		  #--b1 B1               BAM configuration file.
		  #--b2 B2               BAM configuration file.
		  #--s1 S1               FASTQ configuration file.
		  #--s2 S2               FASTQ configuration file.
		  #--od OD               output folder.
		  #-t {paired,single}    readtype, single or paired.
		  #--libType {fr-unstranded,fr-firststrand,fr-secondstrand}
		  #                      Library type. Default is unstranded (fr-unstranded).
		  #                      Use fr-firststrand or fr-secondstrand for strand-
		  #                      specific data.
		  #--readLength READLENGTH
		  #                      The length of each read.
		  #--anchorLength ANCHORLENGTH
		  #                      The anchor length. (default is 1.)
		  #--tophatAnchor TOPHATANCHOR
		  #                      The "anchor length" or "overhang length" used in the
		  #                      aligner. At least “anchor length” NT must be
		  #                      mapped to each end of a given junction. The default is
		  #                      6. (This parameter applies only if using fastq).
		  #--bi BINDEX           The folder name of the STAR binary indexes (i.e., the
		  #                      name of the folder that contains SA file). For
		  #                      example, use ~/STARindex/hg19 for hg19. (Only if using
		  #                      fastq)
		  #--nthread NTHREAD     The number of thread. The optimal number of thread
		  #                      should be equal to the number of CPU core.
		  #--tstat TSTAT         the number of thread for statistical model.
		  #--cstat CSTAT         The cutoff splicing difference. The cutoff used in the
		  #                      null hypothesis test for differential splicing. The
		  #                      default is 0.0001 for 0.01% difference. Valid: 0 ≤
		  #                      cutoff < 1.
		  #--statoff             Turn statistical analysis off.

python /rMATS.4.0.2/rMATS-turbo-Linux-UCS4  --gtf GRCh38_star/${gtf} --b1 `python /rMATS.4.0.2/rMATS-turbo-Linux-UCS4/manifest_to_csl.py manifest.list 1 .` --b2 `python /rMATS.4.0.2/rMATS-turbo-Linux-UCS4/manifest_to_csl.py manifest.list 2 .` --od ${comparison}

aws s3 sync ${comparison} s3://panorama-clk-repro/${project}/${comparison}

