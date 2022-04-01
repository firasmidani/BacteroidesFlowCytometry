#!/usr/bin/env python

# Firas Said Midani
# Start date: 2018-07-10
# Final date: 2022-03-07

# DESCRIPTION pipeline for processing, sampling, and probabilistically summarizing  
#             flow cytometry data for the Bacteroides study (Midani et al. 2022)

###################################
## IMPORT OFF-THE-SHELF PACKAGES ##
###################################

import os
import sys
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

print
print 'PAKCAGES LOADED\n'
print 'pandas\t\t%s' % pd.__version__
print 'matplotlib\t\t%s' %mpl.__version__ 
print 

##############################
## IMPORT IN-HOUSE PACKAGES ##
##############################

ParentPath = "/data/davidlab/users/fsm/bacteroides_20220307/"
PackagePath = ParentPath + 'code/packages/'

sys.path.append(PackagePath)

if not os.path.isdir(PackagePath):  
    sys.exit('Did not find FlowCytometryLibrary.')
else:
    print('Importing FlowCytometryLibrary from ' + PackagePath)

from FlowCytometryLibrary.process import getFormattedTime, readFCS, dataFromFCS, sampleData, dirtyGMM

##########################
## INITIALIZE VARIABLES ##
##########################

print 'Current Time is %s' % getFormattedTime()
print 

####################################
## INITIALIZE FOLDERS FOR RESULTS ##
####################################

ls_dir = ['counts','samples', 'models', 'models_figures','counts/subsets/']
ls_dir = [ParentPath + '/' + ii for ii in ls_dir]

for dir_path in ls_dir: 
    if not os.path.isdir(dir_path): 
        os.mkdir(dir_path)

#######################
## READ MAPPING FILE ##
#######################

mapping = pd.read_csv(ParentPath + '/mapping/mapping.txt',sep='\t',header=0,index_col=0)

print 'Mapping has %s row and %s columns' % (mapping.shape[0],mapping.shape[1])
print 

# Which files to run based on sys.argv[1]
#     Due to the large number of Flow Cytometry files, this code was tweaked to enable
#     parallelization. Samples or files are split into chunks. Below, 'n' or 'numPerChun' 
#     hard-codes the chunk size of 100. The first command-line argument dictates which 
#     chunk size should be analyzed by each execution of this script. Next few lines 
#     therefore identify the exact index of samples (in mapping file) that should be analyzed. 

samples_list = range(mapping.shape[0]);
n = numPerChunk = 100
chunks = [samples_list[i * n:(i+1) * n] for i in range((len(samples_list) + n - 1) // n)];
toRun = chunks[int(sys.argv[1])]

###########################
## INITIALIZE PARAMETERS ##
###########################

# Initialize data frame for storing sample counts
count_df = pd.DataFrame(index=mapping.index,columns=['EventsCount'])

# Whether to analyze flow cytometry events using "A" (area of signals) or "H" (height of signals).
dtype = 'H'

# Which flow cytometry variables to include in analysis
columns = ['{}-{}'.format(cc,dtype) for cc in ['SYTO','RFP','GFP','SSC','FSC']]

# Number of maximum clusters based on sample type
max_nc_dict = {'BO':2,'BF':2,'BT':2,'BV':4,'MIX':10,'BMM':2,'PBS':0}

print "The following mapping rows (indices) will be analyzed"
print toRun
print 

sys.stdout.flush()

####################
## READ FLOW DATA ##
####################

for idx,row in mapping.iloc[toRun,:].iterrows():

    # define absolute file name and path
    FilePath = '%s/data/%s/%s/%s' % (ParentPath,row['Folder'],row['SubFolder'],idx); 
    FileName = idx;
    FileOutput_Sample = '%s/samples/%s.txt' % (ParentPath,idx);

    # describe sample
    sample_info = "_".join([str(ii) for ii in row.loc[['Substrate','Species','TimePoint']].values])

    # communicate
    print 'Processing\t%s\t%s' % (FileName,sample_info)
    sys.stdout.flush()

    # read file
    FCS = readFCS(FilePath);
    data = dataFromFCS(FCS,ZeroFloor=True);

    # sample data 
    data, count = sampleData(data,N=10000,SYTO=450,sample=True);
    print(count)

    # record count of events in each sample
    count_df.loc[idx,'EventsCount'] = count;

    # save reduced version of FCS
    data = data.loc[:,columns]
    data.to_csv(FileOutput_Sample,sep='\t',header=True,index=True);

    # select maximum number of components
    max_nc = max_nc_dict[row['Species']];

    # run GMM, then summarize and visualize
    summary, fig, ax = dirtyGMM(data,max_nc=max_nc,dtype=dtype)

    # save GMM inference predictions
    summary.to_csv('%s/models/%s.txt' % (ParentPath,FileName), sep='\t',header=True,index=True);

    # save GMM visualization
    fig
    fig.subplots_adjust(left=0.25,bottom=0.15,top=0.9)
    fig.savefig('%s/models_figures/%s.pdf' % (ParentPath,FileName), filetype='pdf', bbox_inches='tight')
    plt.close(fig)
#endfor

# save pandas.DataFrame of events count
count_path = '%s/counts/subsets/EventsCount-InferComponents-%s.txt' % (ParentPath,sys.argv[1])
count_df.dropna().to_csv(count_path,sep='\t',header=True,index=True)

print
print '~~~~~ done ~~~~~'

sys.stdout.flush()
