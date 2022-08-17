from ctypes import sizeof
from operator import indexOf
from unicodedata import name
import numpy as np
from numpy.core.fromnumeric import size
import pyhf
import uproot
import hist
from hist import Hist
import json
from optparse import OptionParser
try:
    from HiggsAnalysis.CombinedLimit.Datacard import Datacard
except:
    print("Either the docker container has not been created properly or Combine commands have not been mounted. Please fix this and try again")


def addChannels():  ##add each channel to Datacard object, add observed counts
    DC.bins = [channel["name"] for i, channel in enumerate(spec["channels"])]
    for idxc, channel in enumerate(spec["channels"]):
        if channel_bins[channel["name"]] != 1:
            DC.hasShapes = True
    for idxc, channel in enumerate(spec["channels"]):
        DC.exp.update({channel["name"]: {}})
        if DC.hasShapes == False:  ##single bin
            DC.obs.update({channel["name"]: spec["observations"][idxc]["data"][0]})

        else:
            data = sum(spec["observations"][idxc]["data"])
            DC.obs.update({channel["name"]: data})
            h_data = hist.Hist.new.Regular(
                channel_bins[channel["name"]], 0, channel_bins[channel["name"]]
            ).Weight()
            h_data[...] = np.stack(
                [
                    spec["observations"][idxc]["data"],
                    [0 for i in range(channel_bins[channel["name"]])],
                ],
                axis=-1,
            )

            file[channel["name"] + "/data_obs"] = h_data


def addSamples():  ##add sample names and expected counts to the Datacard
    DC.processes = samples

    for idxc, channel in enumerate(spec["channels"]):
        if DC.hasShapes == False:  ##counting
            for idxs, sample in enumerate(channel["samples"]):
                DC.exp[channel["name"]].update({sample["name"]: sample["data"][0]})
        else:  ##shapes
            DC.shapeMap.update({channel["name"]: {}})
            DC.shapeMap[channel["name"]].update(
                {"data_obs": [options.shapefile, channel["name"] + "/" + "data_obs"]}
            )
            for idxs, sample in enumerate(channel["samples"]):
                data = sum(sample["data"])
                DC.exp[channel["name"]].update({sample["name"]: data})

                mods = [
                    spec["channels"][idxc]["samples"][idxs]["modifiers"][i]["type"]
                    for i, mod in enumerate(
                        spec["channels"][idxc]["samples"][idxs]["modifiers"]
                    )
                ]
                if "histosys" in mods:  ##if shape uncertainty histogram is available
                    DC.shapeMap[channel["name"]].update(
                        {
                            spec["channels"][idxc]["samples"][idxs]["name"]: [
                                options.shapefile,
                                channel["name"]
                                + "/"
                                + spec["channels"][idxc]["samples"][idxs]["name"],
                                channel["name"]
                                + "/"
                                + spec["channels"][idxc]["samples"][idxs]["name"]
                                + "_$SYSTEMATIC",
                            ]
                        }
                    )
                else:
                    DC.shapeMap[channel["name"]].update(
                        {
                            spec["channels"][idxc]["samples"][idxs]["name"]: [
                                options.shapefile,
                                channel["name"]
                                + "/"
                                + spec["channels"][idxc]["samples"][idxs]["name"],
                            ]
                        }
                    )
                if (
                    "staterror" in mods
                ):  ##if staterror is present, provide variances for histograms
                    h_data = hist.Hist.new.Regular(
                        channel_bins[channel["name"]], 0, channel_bins[channel["name"]]
                    ).Weight()

                    h_data[...] = np.stack(
                        [
                            spec["channels"][idxc]["samples"][idxs]["data"],
                            [
                                i**2
                                for i in spec["channels"][idxc]["samples"][idxs][
                                    "modifiers"
                                ][mods.index("staterror")]["data"]
                            ],
                        ],
                        axis=-1,
                    )
                    DC.binParFlags.update({channel["name"]: True})
                elif "shapesys" in mods:  ##for shapesys, also provide variances
                    h_data = hist.Hist.new.Regular(
                        channel_bins[channel["name"]], 0, channel_bins[channel["name"]]
                    ).Weight()
                    h_data[...] = np.stack(
                        [
                            spec["channels"][idxc]["samples"][idxs]["data"],
                            [
                                i**2
                                for i in spec["channels"][idxc]["samples"][idxs][
                                    "modifiers"
                                ][mods.index("shapesys")]["data"]
                            ],
                        ],
                        axis=-1,
                    )
                    DC.binParFlags.update({channel["name"]: True})
                else:  ##no staterror or shapesys->0 variance
                    h_data = hist.Hist.new.Regular(
                        channel_bins[channel["name"]], 0, channel_bins[channel["name"]]
                    ).Weight()
                    h_data[...] = np.stack(
                        [
                            spec["channels"][idxc]["samples"][idxs]["data"],
                            [0 for i in range(channel_bins[channel["name"]])],
                        ],
                        axis=-1,
                    )
                file[
                    channel["name"]
                    + "/"
                    + spec["channels"][idxc]["samples"][idxs]["name"]
                ] = h_data


def addMods():  ##add systematics to DC
    for im, modifier in enumerate(systs):
        if "normsys" in modifier[1]:  ##normsys
            DC.systs.append(
                (modifier[0], False, "shape?", [], {})
            )  ##write normsys as 'shape?' so that Combine doesn't try to combine normsys and histosys mods of the same name
            for idxc, channel in enumerate(spec["channels"]):
                DC.systs[im][4].update({channel["name"]: {}})
                for idxs, sample in enumerate(channel["samples"]):
                    for i in DC.systs:
                        i[4][channel["name"]].update({sample["name"]: 0.0})
        if "lumi" in modifier[1]:  ##lumi
            DC.systs.append(
                (modifier[0], False, "lnN", [], {})
            )  ##write lumi as lnN since they act the same way on the model
            for idxc, channel in enumerate(spec["channels"]):
                DC.systs[im][4].update({channel["name"]: {}})
                for idxs, sample in enumerate(channel["samples"]):
                    for i in DC.systs:
                        i[4][channel["name"]].update({sample["name"]: 0.0})
        if "histosys" in modifier[1]:  ##histosys
            DC.systs.append((modifier[0], False, "shape", [], {}))
            for idxc, channel in enumerate(spec["channels"]):
                DC.systs[im][4].update({channel["name"]: {}})
                for idxs, sample in enumerate(channel["samples"]):
                    for i in DC.systs:
                        i[4][channel["name"]].update({sample["name"]: 0.0})

    for idxc, channel in enumerate(spec["channels"]):
        for idxs, sample in enumerate(channel["samples"]):
            mods = sample["modifiers"]
            names = []
            for mod in mods:
                names.append(mod["name"])
            for syst in DC.systs:
                name = syst[0]
                type = syst[2]
                if name in names:  ##if systematic name is a modifier for this sample
                    if "shape?" in type:  ##normsys
                        for mod in mods:
                            if mod["type"] == "normsys" and mod["name"] == name:
                                if mod["data"]["lo"] == 0:
                                    syst[4][channel["name"]].update(
                                        {
                                            sample["name"]: str(
                                                mod["data"]["lo"] + 1e-9
                                            )
                                            + "/"
                                            + str(mod["data"]["hi"])
                                        }  ##asymmetric lnN
                                    )
                                elif mod["data"]["hi"] == 0:
                                    syst[4][channel["name"]].update(
                                        {
                                            sample["name"]: str(mod["data"]["lo"])
                                            + "/"
                                            + str(mod["data"]["hi"] + 1e-9)
                                        }  ##asymmetric lnN
                                    )
                                else:
                                    syst[4][channel["name"]].update(
                                        {
                                            sample["name"]: str(mod["data"]["lo"])
                                            + "/"
                                            + str(mod["data"]["hi"])
                                        }  ##asymmetric lnN
                                    )
                    if "lnN" in type:  ##lumi only
                        for mod in mods:
                            if mod["type"] == "lumi" and mod["name"] == name:
                                for im, measurement in enumerate(spec["measurements"]):
                                    for ip, param in enumerate(
                                        measurement["config"]["parameters"]
                                    ):
                                        if mod["name"] == param["name"]:
                                            syst[4][channel["name"]].update(
                                                {
                                                    sample["name"]: str(
                                                        param["auxdata"][0]
                                                        - param["sigmas"][0]
                                                    )
                                                    + "/"
                                                    + str(
                                                        param["auxdata"][0]
                                                        + param["sigmas"][0]
                                                    )
                                                }  ##asymmetric lnN
                                            )

                    if "shape" in type:  ##histosys
                        for mod in mods:
                            if mod["type"] == "histosys" and mod["name"] == name:
                                syst[4][channel["name"]].update({sample["name"]: 1.0})
                                hi_data = hist.Hist.new.Regular(
                                    channel_bins[channel["name"]],
                                    0,
                                    channel_bins[channel["name"]],
                                ).Weight()
                                hi_data[...] = np.stack(
                                    [
                                        mod["data"]["hi_data"],
                                        [
                                            0
                                            for i in range(
                                                channel_bins[channel["name"]]
                                            )
                                        ],
                                    ],
                                    axis=-1,
                                )
                                lo_data = hist.Hist.new.Regular(
                                    channel_bins[channel["name"]],
                                    0,
                                    channel_bins[channel["name"]],
                                ).Weight()
                                lo_data[...] = np.stack(
                                    [
                                        mod["data"]["lo_data"],
                                        [
                                            0
                                            for i in range(
                                                channel_bins[channel["name"]]
                                            )
                                        ],
                                    ],
                                    axis=-1,
                                )
                                file[
                                    channel["name"]
                                    + "/"
                                    + spec["channels"][idxc]["samples"][idxs]["name"]
                                    + "_"
                                    + name
                                    + "Up"
                                ] = hi_data
                                file[
                                    channel["name"]
                                    + "/"
                                    + spec["channels"][idxc]["samples"][idxs]["name"]
                                    + "_"
                                    + name
                                    + "Down"
                                ] = lo_data


def addSignal():  ##determine which samples are signal
    measurements = []
    for im, measurement in enumerate(spec["measurements"]):
        measurements.append(measurement["config"]["poi"])
    signalMods = []
    for modifier in modifiers:
        if modifier[0] in measurements:
            signalMods.append(modifier[0])
    for idxc, channel in enumerate(channels):
        for idxs, sample in enumerate(spec["channels"][idxc]["samples"]):
            for i, mod in enumerate(
                spec["channels"][idxc]["samples"][idxs]["modifiers"]
            ):
                for sig in signalMods:
                    if sig == mod["name"]:
                        DC.isSignal.update({sample["name"]: True})
                        DC.signals.append(sample["name"])


def addRateParams():  ##add normfactor mods as rateParams (excluding signal strength)
    measurements = []
    for im, measurement in enumerate(spec["measurements"]):
        measurements.append(measurement["config"]["poi"])
    signalMods = []
    for modifier in modifiers:
        if modifier[0] in measurements:
            signalMods.append(modifier[0])
    print(signalMods)

    for idxc, channel in enumerate(channels):
        for idxs, sample in enumerate((spec["channels"][idxc]["samples"])):
            isSig = False
            for i, mod in enumerate(sample["modifiers"]):
                if mod["name"] in signalMods:
                    isSig = True
            print(sample["name"])
            print(isSig)
            if isSig == False:
                for i, mod in enumerate(
                    spec["channels"][idxc]["samples"][idxs]["modifiers"]
                ):  ##normfactor or shapefactor
                    if "normfactor" in mod["type"] or "shapefactor" in mod["type"]:
                        for im, measurement in enumerate(spec["measurements"]):
                            for ip, param in enumerate(
                                measurement["config"]["parameters"]
                            ):
                                if mod["name"] == param["name"]:
                                    DC.rateParams.update(
                                        {channel + "AND" + sample["name"]: []}
                                    )
                                    DC.rateParams[
                                        channel + "AND" + sample["name"]
                                    ].append([[mod["name"], 1, 0, param["bounds"]], ""])
                                else:
                                    DC.rateParams.update(
                                        {channel + "AND" + sample["name"]: []}
                                    )
                                    DC.rateParams[
                                        channel + "AND" + sample["name"]
                                    ].append([[mod["name"], 1, 0], ""])


def writeDataCard(path):  ##manually write each line of the datacard from DC object
    with open(path, "w") as f:
        f.write("imax " + str(size(DC.bins)) + "\n")
        f.write("jmax " + str(size(DC.processes) - size(DC.isSignal.keys())) + "\n")
        f.write("kmax " + str(size(DC.systs, 0)) + "\n")
        if DC.hasShapes:
            for channel in DC.shapeMap.keys():
                for sample in DC.shapeMap[channel].keys():
                    f.write(
                        "shapes "
                        + sample
                        + "  "
                        + channel
                        + "  "
                        + DC.shapeMap[channel][sample][0]
                        + "  "
                        + DC.shapeMap[channel][sample][1]
                    )
                    if size(DC.shapeMap[channel][sample]) > 2:
                        f.write("  " + DC.shapeMap[channel][sample][2] + "\n")
                    else:
                        f.write("\n")

        f.write("\n---------------------------------\n")
        f.write("bin ")
        for bin in DC.obs.keys():
            f.write(bin + " ")
        f.write("\n")
        f.write("observation ")
        for channel in DC.obs.keys():
            f.write(str(DC.obs[channel]) + " ")
        f.write("\n---------------------------------\n")
        f.write("bin     ")
        for channel in DC.obs.keys():
            for sample in DC.exp[channel].keys():
                f.write(channel + "    ")
        f.write("\n")
        f.write("process     ")
        for channel in DC.bins:
            for sample in DC.exp[channel].keys():
                f.write(sample + "    ")
        f.write("\n")
        f.write("process     ")
        for channel in DC.bins:
            for sample in DC.exp[channel].keys():
                if sample in DC.signals:
                    f.write(str(-1 * DC.processes.index(sample)) + "     ")
                else:
                    f.write(str(DC.processes.index(sample) + 1) + "     ")
        f.write("\n")
        f.write("rate     ")
        for channel in DC.bins:
            for sample in DC.exp[channel].keys():

                f.write(str(DC.exp[channel][sample]) + "     ")
        f.write("\n---------------------------------\n")
        for syst in DC.systs:
            f.write(syst[0] + "  " + syst[2] + "  ")
            for bin in syst[4].keys():
                for sample in DC.exp[bin].keys():
                    if syst[4][bin][sample] != 0:
                        f.write(str(syst[4][bin][sample]) + "  ")
                    else:
                        f.write("-       ")

            f.write("\n")
        f.write("\n---------------------------------\n")
        for cAp in DC.rateParams.keys():
            dir = cAp.split("AND")
            for i in range(size(DC.rateParams[cAp], 0)):
                if size(DC.rateParams[cAp][i][0]) > 3:
                    f.write(
                        str(DC.rateParams[cAp][i][0][0])
                        + " "
                        + "rateParam "
                        + dir[0]
                        + " "
                        + dir[1]
                        + " "
                        + str(DC.rateParams[cAp][i][0][1])
                        + " "
                        + DC.rateParams[cAp][i][0][3]
                    )
                else:
                    f.write(
                        str(DC.rateParams[cAp][i][0][0])
                        + " "
                        + "rateParam "
                        + dir[0]
                        + " "
                        + dir[1]
                        + " "
                        + str(DC.rateParams[cAp][i][0][1])
                    )
                f.write("\n")
        f.write("\n---------------------------------\n")
        for idxc, channel in enumerate(channels):
            if channel in DC.binParFlags.keys():
                if DC.binParFlags[channel] == True:  ##double check to be safe
                    shapesys = False
                    staterror = False
                    for idxs, sample in enumerate(spec["channels"][idxc]["samples"]):
                        if "shapesys" in [
                            mod["type"] for i, mod in enumerate(sample["modifiers"])
                        ]:
                            shapesys = True
                        elif "staterror" in [
                            mod["type"] for i, mod in enumerate(sample["modifiers"])
                        ]:
                            staterror = True
                    if shapesys == True:
                        f.write(
                            channel
                            + " autoMCStats "
                            + str(100000)
                            + " "
                            + str(0)
                            + " "
                            + str(2)
                            + "\n"
                        )
                    if staterror == True:
                        f.write(
                            channel
                            + " autoMCStats "
                            + str(0)
                            + " "
                            + str(0)
                            + " "
                            + str(2)
                            + "\n"
                        )
        f.close()

def main():
    parser = OptionParser()  # add command line args
    parser.add_option(
        "-O", "--out-datacard", dest="outdatacard", default="converted_datacard.txt", help = "desired name of datacard file"
    )
    parser.add_option("-s", "--shape-file", dest="shapefile", default="shapes.root", help = "desired name of shapes file")
    options, args = parser.parse_args()

    with open(args[0]) as serialized:
        spec = json.load(serialized)
    workspace = pyhf.Workspace(spec)
    model = workspace.model()


    channels = model.config.channels
    channel_bins = model.config.channel_nbins
    samples = model.config.samples
    modifiers = model.config.modifiers
    lumi = model.config.parameters

    # generate list of mods without normfactors, shapesys, staterror
    systs = []
    for mod in modifiers:
        if "normfactor" not in mod:
            if "shapesys" not in mod:
                if "staterror" not in mod:
                    systs.append(mod)


    DC = Datacard()


    file = uproot.recreate(options.shapefile)
    addChannels()
    addSamples()
    addMods()
    addSignal()
    addRateParams()
    file.close()
    writeDataCard(options.outdatacard)

if __name__ == "__main__":
    main()
