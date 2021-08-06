import pandas
pandas.set_option('display.max_colwidth',1000)
pandas.set_option('display.max_rows', None)
st = pandas.read_csv("metadata/metadata.csv",
                     dtype={'Run,': object,
                            'ReleaseDate,': object, 'LoadDate,': object, 'spots,': object, 'bases,': object,
                            'spots_with_mates,': object, 'avgLength,': object, 'size_MB,': object, 'AssemblyName,': object,
                            'download_path,': object, 'Experiment,': object, 'LibraryName,': object, 'LibraryStrategy,': object,
                            'LibrarySelection,': object, 'LibrarySource,': object, 'LibraryLayout,': object, 'InsertSize,': object,
                            'InsertDev,': object, 'Platform,': object, 'Model,': object, 'SRAStudy,': object, 'BioProject,': object,
                            'Study_Pubmed_id': object, 'ProjectID,': object, 'Sample,': object, 'BioSample,': object,
                            'SampleType,': object, 'TaxID,': object, 'ScientificName,': object, 'SampleName,': object,
                            'g1k_pop_code,': object, 'source,': object, 'g1k_analysis_group,': object, 'Subject_ID,': object,
                            'Sex,': object, 'Disease,': object, 'Tumor,': object, 'Affection_Status,': object, 'Analyte_Type,': object,
                            'Histological_Type,': object, 'Body_Site,': object, 'CenterName,': object, 'Submission,': object,
                            'dbgap_study_accession,': object, 'Consent,': object, 'RunHash,': object, 'ReadHash': object})

fl = pandas.read_csv("metadata/filenames.txt",sep="\t")
st = st.merge(fl, how='left', on="Run")
sra = pandas.read_csv("metadata/SRP091981.metadata",sep="\t")
#st = st.merge(sra, how='inner', on="Run")
sra = sra.rename(columns={'run_accession': 'Run'})  
sra['stranded']=sra['experiment_title'].str.contains('(stranded)', regex=False)
st = st.merge(sra,how='left', on="Run")

def illuminaRuns():
    return(st.loc[st['Platform'] == 'ILLUMINA']['Run'].tolist())


def pacbioRuns():
    return(st.loc[st['Platform'] == 'PACBIO_SMRT']['Run'].tolist())


def getType(runName):
    return(getProperty(runName, 'Platform'))


def getProperty(runName, column):
    stseries = st.loc[st['Run'] == runName, column]
    if len(stseries) == 1:
        return(stseries.to_string(index=False))  # pandas is ridiculous
    elif len(stseries) > 1:
        raise ValueError("Not expecting to match multiple run names, just 1")
    else:
        raise ValueError("Can't find that run {0}".format(runName))


def getMemory(runName):
    return(int(getProperty(runName, 'size_MB')))


def getECS(runName, units, program):
    if program == 'minimap':
        mb = 32000
    elif program == 'STAR':
        superIntensive = ['SRR5009539', 'SRR5009525',
                          'SRR5009498', 'SRR5009538', 'SRR5009408', 'SRR5009418']
        size = getMemory(runName)
        if runName in superIntensive:
            mb = 192000
        else:
            mb = max(64000, 10*size)
    elif program == 'IsoModule':
        mb = 192000
    elif program == 'samtoolsindex':
        mb = 16000
    elif program == 'samtoolssubsample':
        mb = 16000
    elif program == 'rmats':
        mb = 16000
    else:
        raise ValueError
    if units == 'bytes':
        return(1048576*mb)
    elif units == 'mb':
        return(mb)

# 2  # number of sample
# 3  # number of replicates in sample #1
# sam1_rep1.bam
# sam1_rep2.bam
# sam1_rep3.bam
# 3  # number of replicates in sample #2
# sam2_rep1.bam
# sam2_rep2.bam
# sam2_rep3.bam


def twoSampleComparisonManifest(samp1, samp2, filename, path_prefix=None, platform='ILLUMINA'):
    text_file = open(filename, "w")
    text_file.write("2\n")  # two-way comparison
    for runName in [samp1, samp2]:
        run = getBamsFromSampleName(runName, path_prefix, include_bai=False, platform=platform)
        text_file.write("{0}\n".format(len(run)))
        for replicate in run:
            text_file.write("{0}\n".format(replicate))

# "panorama-clk-repro/SRP091981/


def getBamsFromSampleName(samp, path_prefix=None, include_bai=True, platform=None):
    if include_bai:
        exts = ['bam', 'bam.bai']
    else:
        exts = ['bam']
    platform_ext = {'ILLUMINA':'Aligned.sortedByCoord.out.md','PACBIO_SMRT':'filtered'}
    runs = getRunsFromSampleName(samp,platform=platform)
    bams = []
    #for platform in ['ILLUMINA','PACBIO_SMRT']:
    if path_prefix:
            bams += ["{0}/{1}.{3}.{2}".format(
                path_prefix, run, ext, platform_ext[getPlatformFromRun(run)]) for run in runs for ext in exts]
    else:
        bams += ["{0}.{2}.{1}".format(
            run, ext, platform_ext[getPlatformFromRun(run)]) for run in runs for ext in exts]
    return(bams)

def getPlatformFromRun(run):
    platform=st.loc[(st['Run']==run)]['Platform'].tolist()[0]
    return(platform)

def getFastqsFromSampleName(samp, path_prefix=None, include_bai=True):
    exts = ['_1.fastq.gz','_2.fastq.gz']
    runs = getRunsFromSampleName(samp)
    if path_prefix:
        fastqs = ["{0}/{1}{2}".format(
            path_prefix, replicate, ext) for replicate in runs for ext in exts]
    else:
        fastqs = ["{0}{1}".format(
            replicate, ext) for replicate in runs for ext in exts]
    return(fastqs)

def getRunsFromSampleName(samp,platform=None,stranded=None):
    if stranded is not None:
        #accept either dosage nicknames or the actual sample name
        if platform is not None:
            if samp.lower() in dosageTable:
                return(st.loc[(st['SampleName'].isin(dosageTable[samp.lower()])) & (st['Platform'] == platform) & (st['stranded'] == stranded)]['Run'].tolist())
            else:
                return(st.loc[(st['SampleName']==samp) & (st['Platform'] == platform) & (st['stranded'] == stranded)]['Run'].tolist())
        else:
            if samp.lower() in dosageTable:
                return(st.loc[(st['SampleName'].isin(dosageTable[samp.lower()])) & (st['stranded'] == stranded)]['Run'].tolist())
            else:
                return(st.loc[(st['SampleName']==samp) & (st['stranded'] == stranded)]['Run'].tolist())
    else:
        if platform is not None:
            if samp.lower() in dosageTable:
                return(st.loc[(st['SampleName'].isin(dosageTable[samp.lower()])) & (st['Platform'] == platform) ]['Run'].tolist())
            else:
                return(st.loc[(st['SampleName']==samp) & (st['Platform'] == platform) ]['Run'].tolist())
        else:
            if samp.lower() in dosageTable:
                return(st.loc[(st['SampleName'].isin(dosageTable[samp.lower()])) ]['Run'].tolist())
            else:
                return(st.loc[(st['SampleName']==samp) ]['Run'].tolist())

def getfulldosagename(nickname):
    return(dosageTable[nickname])

def printDoseMintieSymlinks(treatment):
    files=getFastqsFromSampleName(treatment)
    print("mkdir -p mintiesymlinks/{}".format(treatment))
    for f in files:
        print("ln -s ../../SRP091981/{0} mintiesymlinks/{1}/{0}".format(f,treatment))

dosageTable = {'untreated': ['Untreated HCT116'],  '0.5': ['0.5 uM T3 treated HCT116'], '0.1': ['0.1 uM T3 treated HCT116'], '0.05': ['0.05 uM T3 treated HCT116'], '1.0': ['1.0 uM T3 treated HCT116'],'5.0':['5 uM T3 treated HCT116'],
               'untreated184': ['Untreated 184-hTert'], '0.5-184': ['0.5 uM T3 treated 184-hTert'], '1.0-184': ['1.0 uM T3 treated 184-hTert'], '5.0-184': ['5.0 uM T3 treated 184-hTert']}
dosageTable['treated']=dosageTable['0.5']+dosageTable['0.1']+dosageTable['0.05']+dosageTable['1.0']
dosageTable['treated184']=dosageTable['0.5-184']+dosageTable['1.0-184']+dosageTable['5.0-184']

def printAllMintieSymlinks():
    for dose in dosageTable:
        printDoseMintieSymlinks(dose)


def getSUPPAconfig(dose):
    runs = getRunsFromSampleName(dose)
    for run in runs:
        print('  - "{0}"'.format(run))

#metautils.printSuppaComparison("treated","untreated")
#metautils.printSuppaComparison("0.5","untreated")
def printSuppaComparison(experiment,control):
    print("EXPERIMENT:")
    getSUPPAconfig(experiment)
    print("CONTROL:")
    getSUPPAconfig(control)