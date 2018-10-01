import pandas

st = pandas.read_csv("metadata.csv",
dtype={'Run,':object,
'ReleaseDate,':object,'LoadDate,':object,'spots,':object,'bases,':object,
'spots_with_mates,':object,'avgLength,':object,'size_MB,':object,'AssemblyName,':object,
'download_path,':object,'Experiment,':object,'LibraryName,':object,'LibraryStrategy,':object,
'LibrarySelection,':object,'LibrarySource,':object,'LibraryLayout,':object,'InsertSize,':object,
'InsertDev,':object,'Platform,':object,'Model,':object,'SRAStudy,':object,'BioProject,':object,
'Study_Pubmed_id':object,'ProjectID,':object,'Sample,':object,'BioSample,':object,
'SampleType,':object,'TaxID,':object,'ScientificName,':object,'SampleName,':object,
'g1k_pop_code,':object,'source,':object,'g1k_analysis_group,':object,'Subject_ID,':object,
'Sex,':object,'Disease,':object,'Tumor,':object,'Affection_Status,':object,'Analyte_Type,':object,
'Histological_Type,':object,'Body_Site,':object,'CenterName,':object,'Submission,':object,
'dbgap_study_accession,':object,'Consent,':object,'RunHash,':object,'ReadHash':object})

def illuminaRuns():
    return(st.loc[st['Platform']=='ILLUMINA']['Run'].tolist())

def pacbioRuns():
    return(st.loc[st['Platform']=='PACBIO_SMRT']['Run'].tolist())
    
def getType(runName):
    getProperty(runName,'Platform')
    
def getProperty(runName, column):
    stseries=st.loc[st['Run'] == runName, column]
    if len(stseries)==1:
        return(stseries.to_string(index=False)) #pandas is ridiculous
    elif len(stseries)>1:
        raise ValueError("Not expecting to match multiple run names, just 1")
    else:
        raise ValueError("Can't find that run {0}".format(runName))

def getMemory(runName):
    return(int(getProperty(runName,'size_MB')))


def getECS(runName,units,program):
    if program=='minimap':
        mb = 32000
    elif program=='STAR':
        superIntensive = ['SRR5009539', 'SRR5009525', 'SRR5009498', 'SRR5009538', 'SRR5009408', 'SRR5009418']
        size = getMemory(runName)
        if runName in superIntensive:
            mb = 192000
        else:
            mb = max(64000,10*size)
    else:
        raise ValueError
    if units=='bytes':
        return(1048576*mb)
    elif units=='mb':
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
def twoSampleComparisonManifest(samp1,samp2,filename):
    run1 = getBamsFromSampleName(samp1)
    run2 = getBamsFromSampleName(samp2)
    text_file = open(filename, "w")
    text_file.write("2\n") #two-way comparison
    for runName in [samp1,samp2]:
        run = st.loc[st['SampleName']==runName]['Run'].tolist()
        text_file.write("{0}\n".format(len(run)))
        for replicate in run:
            text_file.write("{0}\n".format(replicate))

def getBamsFromSampleName(samp):
    run = getRunsFromSampleName(samp)
    bams = ["{0}.Aligned.sortedByCoord.out.bam".format(replicate) for replicate in run]
    return(bams)

def getRunsFromSampleName(samp):
    return(st.loc[(st['SampleName']==samp) & (st['Platform']=='ILLUMINA')]['Run'].tolist())