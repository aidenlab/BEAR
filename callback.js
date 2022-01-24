function launchSpec(dataProvider)
{
    var ret = {
        commandLine: [ "python3", "POLAR-BEAR/basespace_app/app_handler.py"],
        containerImageId: "peradastra/polar_bear_native_app:betaTesting",
        Options: [  "bsfs.enabled=true"]
    };
    return ret;
}