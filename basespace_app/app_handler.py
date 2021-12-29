import os
import fnmatch
import json


# load json file
## define JSON format of the "AppSession.json" for the app session
def metadatajson():
    json_file = json.dumps(
        {
            "Name"          : "",
            "Description"   : "",
            "HrefAppSession": "",
            "Properties"    : [
                {
                    "Type" : "sample[]",
                    "Name" : "Input.Samples",
                    "Items": [

                    ]
                }
            ]
        }
    )
    return json_file

## open and load "AppSession.json" for the app session
jsonfile = open('data/input/AppSession.json')
jsonObject = json.load(jsonfile)

# get info necessary to
# determine the number of properties
numberOfPropertyItems = len(jsonObject['Properties']['Items'])

sampleID = []
projectID = []
R1files = []
R2files = []
for index in range(numberOfPropertyItems):
    if jsonObject['Properties']['Items'][index]['Name'] == 'Input.Projects':
        projectID = jsonObject['Properties']['Items'][index]['Items'][0]['Id']
    if jsonObject['Properties']['Items'][index]['Name'] == 'Input.Samples':
        for sample in range(len(jsonObject['Properties']['Items'][index]['Items'])):
            sampleID.append(jsonObject['Properties']['Items'][index]['Items'][sample]['Id'])
            sampleDir = '/data/input/samples/%s' % (sampleID[sample])
            for root, dirs, files in os.walk(sampleDir[sample]):
                R1files.append(fnmatch.filter(files, '*_R1_*'))
                R2files.append(fnmatch.filter(files, '*_R2_*'))

sampleOutDir = '/data/output/appresults/%s/%s' %(projectID,"output_here")
os.system('mkdir -p "%s"' %(sampleOutDir))

file = '%s/parameters.csv' %(sampleOutDir)
outFile = open(file ,'w')
outFile.write('This is a test.')
outFile.close()
