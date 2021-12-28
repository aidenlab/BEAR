import os
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
    numberOfPropertyItems = len(json['Properties']['Items'])

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
    bearOutDir = '/data/output/appresults/%s/%s' % (projectID, "polar_bear")

    return sampleDir, projectID, bearOutDir

def make_output_folder_in_basespace(bearOutDir):
    # create output folder for BEAR pipeline in the project folder
    os.system('mkdir -p "%s"' %(bearOutDir))

def write_to_parameter_file(bearOutDir, parameterData):
    # create "parameters.csv" file in output folder for BEAR
    file = '%s/parameters.csv' %(bearOutDir)
    parameterFile = open(file ,'w')
    parameterFile.write(parameterData)
    parameterFile.close()

def run_polar_bear():
    # bear_command = ["bash", "POLAR-BEAR/run_polar_bear_eua_pipline.sh", "-d", topDir]
    bear_command = ["bash", "POLAR-BEAR/run_polar_bear_eua_pipline.sh", "-h"]
    bear_out = subprocess.call(bear_command)
    return bear_out


if __name__ == "__main__":
    jsonAppObject = load_appsession_json()
    sampleDir, projectID, bearOutDir = get_app_data_from_json(jsonAppObject)
    make_output_folder_in_basespace(bearOutDir)
    bear_out = run_polar_bear()
    write_to_parameter_file(bearOutDir, bear_out)

