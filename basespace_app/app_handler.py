import os
import subprocess
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

def load_appsession_json(jsonPath='/data/input/AppSession.json'):
    # open and load "AppSession.json" for the app session
    jsonfile = open(jsonPath)
    jsonAppObject = json.load(jsonfile)
    return jsonAppObject

def get_app_data_from_json(jsonObject):
    # get info necessary to
    # determine the number of properties
    numberOfPropertyItems = len(jsonObject['Properties']['Items'])

    # get "projectID" and "sampleID" information
    projectID = []
    sampleID = []
    for index in range(numberOfPropertyItems):
        if jsonObject['Properties']['Items'][index]['Name'] == 'Input.Projects':
            projectID = jsonObject['Properties']['Items'][index]['Items'][0]['Id']
        if jsonObject['Properties']['Items'][index]['Name'] == 'Input.Samples':
            for sample in range(len(jsonObject['Properties']['Items'][index]['Items'])):
                sampleID.append(jsonObject['Properties']['Items'][index]['Items'][sample]['Id'])

    # set sample directory containing FASTQs
    sampleDir = '/data/input/samples/%s' % (sampleID[sample])

    # set output directory
    bearOutDir = '/data/output/appresults/%s/%s' %(projectID,"polar-bear-fda-eua")

    return sampleDir, projectID, bearOutDir


def create_file_trigger_for_bear_pipeline(bearOutDir):
    # create "parameters.csv" file in output folder for BEAR
    parameterFile = open("POLAR-BEAR/native.app.txt" ,'w')
    parameterFile.write(bearOutDir)
    parameterFile.close()

def create_output_folder(bearOutDir):
    # create "parameters.csv" file in output folder for BEAR
    os.system('mkdir -p "%s"' %(bearOutDir))


def run_polar_bear():
    # create & open log file for bear pipeline
    bear_pipeline_log = open('data/logs/polar_bear_eua_pipline_log.txt', 'w')

    # run bear command
    bear_command = ["bash", "POLAR-BEAR/run_polar_bear_eua_pipline.sh", "-d", "POLAR-BEAR/test"]
    bear_out = subprocess.call(bear_command, stdout=bear_pipeline_log)

    # close log file for bear pipeline
    bear_pipeline_log.close()


if __name__ == "__main__":
    jsonAppObject = load_appsession_json()
    sampleDir, projectID, bearOutDir = get_app_data_from_json(jsonAppObject)
    create_file_trigger_for_bear_pipeline(bearOutDir)
    create_output_folder(bearOutDir)
    run_polar_bear()


#


