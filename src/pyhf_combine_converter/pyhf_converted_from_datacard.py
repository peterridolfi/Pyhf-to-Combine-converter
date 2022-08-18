import string
from numpy.core.fromnumeric import size
import uproot
import json
from optparse import OptionParser

try:
    import HiggsAnalysis.CombinedLimit.DatacardParser as DP
    from HiggsAnalysis.CombinedLimit.Datacard import Datacard
except:
    print(
        "Either the docker container has not been created properly or Combine commands have not been mounted. Please fix this and try again."
    )


def getShapeFile(
    shapeMap: dict, channel, sample
) -> string:  ##get shape file of specific sample
    file = ""
    if not channel in shapeMap.keys():
        if "*" in shapeMap.keys():
            if not sample in shapeMap["*"]:
                if "*" in shapeMap["*"].keys():
                    file = shapeMap["*"]["*"][0]
            else:
                file = shapeMap["*"][sample][0]
    elif not sample in shapeMap[channel]:
        if "*" in shapeMap[channel].keys():
            file = shapeMap[channel]["*"][0]
    else:
        file = shapeMap[channel][sample][0]
    return file


def getHistPath(
    shapeMap: dict, channel, sample
) -> string:  ##get path to sample histogram
    path = ""
    if not channel in shapeMap.keys():
        if "*" in shapeMap.keys():
            if not sample in shapeMap["*"]:
                if "*" in shapeMap["*"].keys():
                    path = (
                        shapeMap["*"]["*"][1]
                        .replace("$CHANNEL", channel)
                        .replace("$PROCESS", sample)
                    )
            else:
                path = shapeMap["*"][sample][1].replace("$CHANNEL", channel)
    elif not sample in shapeMap[channel]:
        if "*" in shapeMap[channel].keys():
            path = shapeMap[channel]["*"][1].replace("$PROCESS", sample)
    else:
        path = shapeMap[channel][sample][1]
    return path


def getUncertPath(
    shapeMap: dict, channel, sample, name
) -> string:  ##get path to shape uncertainty if it exists
    path = ""
    if not channel in shapeMap.keys():
        if "*" in shapeMap.keys():
            if not sample in shapeMap["*"]:
                if "*" in shapeMap["*"].keys():
                    path = (
                        shapeMap["*"]["*"][2]
                        .replace("$CHANNEL", channel)
                        .replace("$PROCESS", sample)
                        .replace("$SYSTEMATIC", name)
                    )
            else:
                path = (
                    shapeMap["*"][sample][2]
                    .replace("$CHANNEL", channel)
                    .replace("SYSTEMATIC", name)
                )
    elif not sample in shapeMap[channel]:
        if "*" in shapeMap[channel].keys():
            path = (
                shapeMap[channel]["*"][2]
                .replace("$PROCESS", sample)
                .replace("$SYSTEMATIC", name)
            )
    else:
        path = shapeMap[channel][sample][2].replace("$SYSTEMATIC", name)
    return path


def getHist(shapeMap: dict, channel, sample):  ##get actual sample histogram
    file = uproot.open(getShapeFile(shapeMap, channel, sample))
    hist = file[getHistPath(shapeMap, channel, sample)]
    return hist


def getUncertUp(
    shapeMap: dict, channel, sample, name
):  ##get shape uncertainty shifted up
    file = uproot.open(getShapeFile(shapeMap, channel, sample))
    hist = file[getUncertPath(shapeMap, channel, sample, name) + "Up"]
    return hist


def getUncertDown(
    shapeMap: dict, channel, sample, name
):  ##get shape uncertaint shifted down
    file = uproot.open(getShapeFile(shapeMap, channel, sample))
    hist = file[getUncertPath(shapeMap, channel, sample, name) + "Down"]
    return hist


def addChannels(spec: dict, data_card):
    """
    Add each channel and associated observation to the spec
    """
    if data_card.hasShapes:
        for idxc, channel in enumerate(channels):
            spec["channels"].append({"name": channel, "samples": []})
            hist = getHist(data_card.shapeMap, channel, "data_obs")
            data = hist.values().tolist()
            spec["observations"].append({"name": channel, "data": data})
    else:
        for idxc, channel in enumerate(channels):
            spec["channels"].append({"name": channel, "samples": []})
            spec["observations"].append({"name": channel, "data": [observations[idxc]]})


def addSamples(spec: dict, data_card):
    """
    Add sample names and expected values to spec
    """
    if data_card.hasShapes:
        for idxc, channel in enumerate(channels):
            for idxs, sample in enumerate(samples):
                hist = getHist(data_card.shapeMap, channel, sample)
                data = hist.values().tolist()
                if sample in exp_values[channel].keys():
                    spec["channels"][idxc]["samples"].append(
                        {
                            "name": sample,
                            "data": data,
                            "modifiers": [],
                        }
                    )
    else:
        for idxc, channel in enumerate(channels):
            for idxs, sample in enumerate(samples):
                if sample in exp_values[channel].keys():
                    spec["channels"][idxc]["samples"].append(
                        {
                            "name": sample,
                            "data": [exp_values[channel][sample]],
                            "modifiers": [],
                        }
                    )


def addMeasurements(spec: dict, data_card):
    """
    Add signal measurements to spec
    """
    for idxc, channel in enumerate(channels):
        for idxs, sample in enumerate(samples):
            if sig[sample] == True:
                spec["measurements"].append(
                    {
                        "name": "Measurement_" + sample,
                        "config": {"poi": "mu_" + sample, "parameters": []},
                    }
                )
                if channel + "AND" + sample in data_card.rateParams.keys():
                    for param in data_card.rateParams[channel + "AND" + sample]:
                        spec["measurements"][len(spec["measurements"] - 1)]["config"][
                            "parameters"
                        ].append(
                            {
                                "name": param[0][0],
                                "fixed": False,
                                "inits": [param[0][1]],
                                "bounds": [param[0][3]],
                            }
                        )


def addNormFactor(spec: dict, data_card):
    """
    Add all normfactor modifiers to spec
    """
    for idxc, channel in enumerate(channels):
        for idxs, sample in enumerate(samples):
            if (channel + "AND" + sample) in data_card.rateParams.keys():
                name = data_card.rateParams[channel + "AND" + sample][0][0][0]
                spec["channels"][idxc]["samples"][idxs]["modifiers"].append(
                    {"name": name, "type": "normfactor", "data": None}
                )
            elif sig[sample] == True:  ##add signal strength modifier
                spec["channels"][idxc]["samples"][idxs]["modifiers"].append(
                    {"name": "mu_" + sample, "type": "normfactor", "data": None}
                )


def addMods(spec: dict, data_card):
    """
    Add systematics as modifiers
    """
    for syst in mods:
        name = syst[0]
        mod_type = syst[2]
        if mod_type == "lnN":  ##normsys
            for idxc, channel in enumerate(channels):
                for idxs, sample in enumerate(exp_values[channel].keys()):
                    if sample in exp_values[channel].keys():
                        if syst[4][channel][sample] != 0:
                            if size(syst[4][channel][sample]) == 1:
                                spec["channels"][idxc]["samples"][idxs][
                                    "modifiers"
                                ].append(
                                    {
                                        "name": name,
                                        "type": "normsys",
                                        "data": {
                                            "hi": syst[4][channel][sample],
                                            "lo": 1 / (syst[4][channel][sample]),
                                        },
                                    }
                                )
                            else:
                                data = [i for i in syst[4][channel][sample]]
                                spec["channels"][idxc]["samples"][idxs][
                                    "modifiers"
                                ].append(
                                    {
                                        "name": name,
                                        "type": "normsys",
                                        "data": {
                                            "hi": data[1],
                                            "lo": data[0],
                                        },
                                    }
                                )

        elif "shape" in mod_type:  ##histosys
            for idxc, channel in enumerate(channels):
                for idxs, sample in enumerate(samples):
                    if syst[4][channel][sample] != 0:
                        if sample in exp_values[channel].keys():
                            histUp = getUncertUp(
                                data_card.shapeMap, channel, sample, name
                            )
                            histDown = getUncertDown(
                                data_card.shapeMap, channel, sample, name
                            )
                            hi_data = histUp.values().tolist()
                            lo_data = histDown.values().tolist()
                            data = spec["channels"][idxc]["samples"][idxs]["data"]
                            hi = 0
                            lo = 0
                            nom = 0
                            for i in hi_data:
                                hi = hi + i
                            for i in lo_data:
                                lo = lo + i
                            for i in data:
                                nom = nom + i
                            spec["channels"][idxc]["samples"][idxs][
                                "modifiers"
                            ].append(  ##add histosys modifier for shape
                                {
                                    "name": name,
                                    "type": "histosys",
                                    "data": {"hi_data": hi_data, "lo_data": lo_data},
                                }
                            )
                            spec["channels"][idxc]["samples"][idxs][
                                "modifiers"
                            ].append(  ##add additional normsys to account for shape/norm split
                                {
                                    "name": name,
                                    "type": "normsys",
                                    "data": {"hi": nom / hi, "lo": nom / lo},
                                }
                            )

        else:  ##lnU or gmN
            raise NotImplementedError

    for idxc, channel in enumerate(channels):  ##staterror/shapesys

        if channel in data_card.binParFlags.keys():
            if data_card.binParFlags[channel][0] == 0:  # if BB lite enabled
                for idxs, sample in enumerate(samples):
                    if sample in exp_values[channel].keys():
                        hist = getHist(data_card.shapeMap, channel, sample)
                        err = hist.errors().tolist()
                        spec["channels"][idxc]["samples"][idxs]["modifiers"].append(
                            {"name": "my_stat_err", "type": "staterror", "data": err}
                        )

            elif (
                data_card.binParFlags[channel][0] > 0
            ):  ##if user wants total BB (shapesys)
                for idxs, sample in enumerate(samples):
                    if sample in exp_values[channel].keys():
                        hist = getHist(data_card.shapeMap, channel, sample)
                        err = hist.errors().tolist()
                        spec["channels"][idxc]["samples"][idxs]["modifiers"].append(
                            {
                                "name": "my_shapesys" + channel + sample,
                                "type": "shapesys",
                                "data": err,
                            }
                        )


def main():
    parser = OptionParser()
    DP.addDatacardParserOptions(parser)
    parser.add_option(
        "-O",
        "--out-file",
        dest="outfile",
        default="converted_workspace.json",
        help="desired name of JSON file",
    )
    options, args = parser.parse_args()  # add command line args

    data_card = Datacard()  # create Datacard object
    with open(args[0]) as dc_file:
        data_card = Datacard()
        data_card = DP.parseCard(file=dc_file, options=options)

    channels = [channel for channel in data_card.bins]
    observations = [obs for channel, obs in data_card.obs.items()]
    samples = [sample for sample in data_card.processes]
    exp_values = data_card.exp
    sig = data_card.isSignal
    mods = data_card.systs
    spec = {"channels": [], "observations": [], "measurements": [], "version": "1.0.0"}

    # convert to JSON spec
    addChannels(spec, data_card)
    addSamples(spec, data_card)
    addMeasurements(spec, data_card)
    addNormFactor(spec, data_card)
    addMods(spec, data_card)

    with open(options.outfile, "w") as file:
        file.write(json.dumps(spec, indent=2))


if __name__ == "__main__":
    main()
