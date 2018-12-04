set -euo pipefail
cd $TMPDIR

echo "downloading bams.."
for f in ${untreated}; do aws s3 cp --quiet s3://$f .; done;
for f in ${treated}; do aws s3 cp  --quiet s3://$f .; done;

echo "for f in ${untreated}; do aws s3 cp s3://$f .; done; for f in ${treated}; do aws s3 cp s3://$f .; done;"
aws s3 cp s3://${manifest} manifest.list
echo "aws s3 cp s3://${manifest} manifest.list"

pwd
ls -alt .

echo "downloading gtf.."
mkdir GRCh38_star
aws s3 cp  --quiet s3://panorama-refs/GRCh38_star/${gtf} GRCh38_star/

echo "download rmats results"
aws s3 sync s3://panorama-clk-repro/${project}/${comparison} results
echo "aws s3 sync s3://panorama-clk-repro/${project}/${comparison} ."
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

echo "running rmats2sashimi.."
echo "rmats2sashimiplot --b1 <( python2.7 /manifest_to_csl.py manifest.list 1 . ) \
                  --b2 <( python2.7 /manifest_to_csl.py manifest.list 2 . ) \
                  -t SE \
                  -e SE.MATS.JC.txt \
                  --l1 ${sample1} \
                  --l2 ${sample2} \
                  --exon_s 1 \
                  --intron_s 5 \
                  -o sashimi"

  #-t eventType	                Type of event from rMATS result used in the 											analysis.
  #                                  eventType is 'SE', 'A5SS', 'A3SS', 'MXE' or 'RI'.
  #                                  'SE' is for skipped exon events, 'A5SS' is for
  #                                  alternative 5' splice site events, 'A3SS' is for
  #                                  alternative 3' splice site events, 'MXE' is for
  #                                  mutually exclusive exons events and 'RI' is for
  #                                  retained intron events (Only if using rMATS format
  #                                  result as event file).
  #  -e eventsFile	                The rMATS output event file (Only if using rMATS
  #                                  format result as event file).
                                    
#python2.7 /rMATS.4.0.2/rMATS-turbo-Linux-UCS4/rmats.py  --gtf GRCh38_star/${gtf} --b1 <( python2.7 /rMATS.4.0.2/rMATS-turbo-Linux-UCS4/manifest_to_csl.py manifest.list 1 . ) --b2 <( python2.7 /rMATS.4.0.2/rMATS-turbo-Linux-UCS4/manifest_to_csl.py manifest.list 2 . ) --od ${comparison}

#$rmats2sashimiplot --s1 s1_rep1.sam[,s1_rep2.sam]* --s2 s2.rep1.sam[,s2.rep2.sam]* -t eventType -e eventsFile --l1 SampleLabel1 --l2 SampleLabel2 --exon_s exonScale --intron_s intronScale -o outDir

mkdir -p sashimi/SE sashimi/A5SS sashimi/A3SS sashimi/MXE sashimi/RI

rmats2sashimiplot --b1 `python2.7 /manifest_to_csl.py manifest.list 1 .` \
                  --b2 `python2.7 /manifest_to_csl.py manifest.list 2 .` \
                  -t SE \
                  -e results/SE.MATS.JC.txt \
                  --l1 ${sample1} \
                  --l2 ${sample2} \
                  --exon_s 1 \
                  --intron_s 5 \
                  -o sashimi/SE

rmats2sashimiplot --b1 `python2.7 /manifest_to_csl.py manifest.list 1 .` \
                  --b2 `python2.7 /manifest_to_csl.py manifest.list 2 .` \
                  -t A5SS \
                  -e results/A5SS.MATS.JC.txt \
                  --l1 ${sample1} \
                  --l2 ${sample2} \
                  --exon_s 1 \
                  --intron_s 5 \
                  -o sashimi/A5SS

rmats2sashimiplot --b1 `python2.7 /manifest_to_csl.py manifest.list 1 .` \
                  --b2 `python2.7 /manifest_to_csl.py manifest.list 2 .` \
                  -t A3SS \
                  -e results/A3SS.MATS.JC.txt \
                  --l1 ${sample1} \
                  --l2 ${sample2} \
                  --exon_s 1 \
                  --intron_s 5 \
                  -o sashimi/A3SS

rmats2sashimiplot --b1 `python2.7 /manifest_to_csl.py manifest.list 1 .` \
                  --b2 `python2.7 /manifest_to_csl.py manifest.list 2 .` \
                  -t MXE \
                  -e results/MXE.MATS.JC.txt \
                  --l1 ${sample1} \
                  --l2 ${sample2} \
                  --exon_s 1 \
                  --intron_s 5 \
                  -o sashimi/MXE

rmats2sashimiplot --b1 `python2.7 /manifest_to_csl.py manifest.list 1 .` \
                  --b2 `python2.7 /manifest_to_csl.py manifest.list 2 .` \
                  -t RI \
                  -e results/RI.MATS.JC.txt \
                  --l1 ${sample1} \
                  --l2 ${sample2} \
                  --exon_s 1 \
                  --intron_s 5 \
                  -o sashimi/RI
#rmats2sashimiplot --s1 ./testData/S1.R1.test.sam,./testData/S1.R2.test.sam,./testData/S1.R3.test.sam --s2 ./testData/S2.R1.test.sam,./testData/S2.R2.test.sam,./testData/S2.R3.test.sam 
#-t SE -e ./testData/MATS_output/test_PC3E_GS689.SE.MATS.events.txt --l1 PC3E --l2 GS689 --exon_s 1 --intron_s 5 -o test_events_output  

aws s3 sync sashimi s3://panorama-clk-repro/${project}/${comparison}/sashimi
