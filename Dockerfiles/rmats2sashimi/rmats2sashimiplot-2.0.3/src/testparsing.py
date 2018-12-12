import os
import sys
import argparse
from rmats2sashimiplot.rmats2sashimiplot  import plot_with_eventsfile
class XClass( object ):
   def __init__( self ):
       self.myAttr= None

options = XClass()

options.events_file="SE.MATS.JC.filtered.txt"
options.sashimi_path='sashimi'
options.event_type='SE'
options.out_dir='outdir'
options.b1='/home/leipzig/CLK_reproduction/Dockerfiles/rmats2sashimi/rmats2sashimiplot-2.0.3/src/SRR5009496.Aligned.sortedByCoord.out.bam,/home/leipzig/CLK_reproduction/Dockerfiles/rmats2sashimi/rmats2sashimiplot-2.0.3/src/SRR5009521.Aligned.sortedByCoord.out.bam,/home/leipzig/CLK_reproduction/Dockerfiles/rmats2sashimi/rmats2sashimiplot-2.0.3/src/SRR5009474.Aligned.sortedByCoord.out.bam,/home/leipzig/CLK_reproduction/Dockerfiles/rmats2sashimi/rmats2sashimiplot-2.0.3/src/SRR5009406.Aligned.sortedByCoord.out.bam'
options.b2='/home/leipzig/CLK_reproduction/Dockerfiles/rmats2sashimi/rmats2sashimiplot-2.0.3/src/SRR5009487.Aligned.sortedByCoord.out.bam,/home/leipzig/CLK_reproduction/Dockerfiles/rmats2sashimi/rmats2sashimiplot-2.0.3/src/SRR5009515.Aligned.sortedByCoord.out.bam,/home/leipzig/CLK_reproduction/Dockerfiles/rmats2sashimi/rmats2sashimiplot-2.0.3/src/SRR5009377.Aligned.sortedByCoord.out.bam,/home/leipzig/CLK_reproduction/Dockerfiles/rmats2sashimi/rmats2sashimiplot-2.0.3/src/SRR5009459.Aligned.sortedByCoord.out.bam'
options.exon_s=1
options.intron_s=1
options.font_size=8
options.hide_number=False
options.min_counts=0
options.text_background=True
options.group_info=None
options.color = None
options.l1='untreated'
options.l2='treated'
plot_with_eventsfile(options)