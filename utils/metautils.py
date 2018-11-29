import pandas

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


def twoSampleComparisonManifest(samp1, samp2, filename):
    text_file = open(filename, "w")
    text_file.write("2\n")  # two-way comparison
    for runName in [samp1, samp2]:
        run = getBamsFromSampleName(runName, include_bai=False)
        text_file.write("{0}\n".format(len(run)))
        for replicate in run:
            text_file.write("{0}\n".format(replicate))

# "panorama-clk-repro/SRP091981/


def getBamsFromSampleName(samp, include_s3=None, include_bai=True):
    if include_bai:
        exts = ['bam', 'bam.bai']
    else:
        exts = ['bam']
    runs = getRunsFromSampleName(samp)
    if include_s3:
        bams = ["{0}/{1}.Aligned.sortedByCoord.out.{2}".format(
            include_s3, replicate, ext) for replicate in runs for ext in exts]
    else:
        bams = ["{0}.Aligned.sortedByCoord.out.{1}".format(
            replicate, ext) for replicate in runs for ext in exts]
    return(bams)


def getRunsFromSampleName(samp):
    #accept either dosage nicknames or the actual sample name
    if samp.lower() in dosageTable:
        return(st.loc[(st['SampleName'].isin(dosageTable[samp.lower()])) & (st['Platform'] == 'ILLUMINA')]['Run'].tolist())
    else:
        return(st.loc[(st['SampleName']==samp) & (st['Platform'] == 'ILLUMINA')]['Run'].tolist())


def getfulldosagename(nickname):
    return(dosageTable[nickname])


dosageTable = {'untreated': ['Untreated HCT116'],  '0.5': ['0.5 uM T3 treated HCT116'], '0.1': ['0.1 uM T3 treated HCT116'], '0.05': ['0.05 uM T3 treated HCT116'], '1.0': ['1.0 uM T3 treated HCT116'],'5.0':['5 uM T3 treated HCT116'],
               'untreated184': ['Untreated 184-hTert'], '0.5-184': ['0.5 uM T3 treated 184-hTert'], '1.0-184': ['1.0 uM T3 treated 184-hTert'], '5.0-184': ['5.0 uM T3 treated 184-hTert']}
dosageTable['treated']=dosageTable['0.5']+dosageTable['0.1']+dosageTable['0.05']+dosageTable['1.0']