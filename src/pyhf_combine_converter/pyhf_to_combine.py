import json
from ctypes import sizeof
from operator import indexOf
from optparse import OptionParser
from unicodedata import name

import hist
import numpy as np
import pyhf
import uproot
from hist import Hist
from numpy.core.fromnumeric import size

try:
    from HiggsAnalysis.CombinedLimit.Datacard import Datacard
except:
    print(
        "Either the docker container has not been created properly or Combine commands have not been mounted. Please fix this and try again."
    )

__all__ = ["pyhf_convert_to_datacard"]


def __dir__():
    return __all__


def addChannels(file, spec, data_card, channel_bins):
    """
    Add each channel to Datacard object, add observed counts
    """
    data_card.bins = [channel["name"] for channel in spec["channels"]]
    for channel in spec["channels"]:
        if channel_bins[channel["name"]] != 1:
            data_card.hasShapes = True
    for idxc, channel in enumerate(spec["channels"]):
        data_card.exp.update({channel["name"]: {}})
        # single bin
        if data_card.hasShapes == False:
            data_card.obs.update(
                {channel["name"]: spec["observations"][idxc]["data"][0]}
            )
        else:
            data = sum(spec["observations"][idxc]["data"])
            data_card.obs.update({channel["name"]: data})
            h_data = hist.Hist.new.Regular(
                channel_bins[channel["name"]], 0, channel_bins[channel["name"]]
            ).Weight()

            h_data[...] = np.stack(
                [
                    spec["observations"][idxc]["data"],
                    [0 for _ in range(channel_bins[channel["name"]])],
                ],
                axis=-1,
            )

            file[channel["name"] + "/data_obs"] = h_data


def addSamples(file, spec, data_card, channel_bins, samples, shapefile):
    """
    Add sample names and expected counts to the Datacard
    """
    data_card.processes = samples

    for idxc, channel in enumerate(spec["channels"]):
        if data_card.hasShapes == False:  # counting
            for sample in channel["samples"]:
                data_card.exp[channel["name"]].update(
                    {sample["name"]: sample["data"][0]}
                )
        else:
            # shapes
            data_card.shapeMap.update({channel["name"]: {}})
            data_card.shapeMap[channel["name"]].update(
                {"data_obs": [shapefile, channel["name"] + "/" + "data_obs"]}
            )
            for idxs, sample in enumerate(channel["samples"]):
                data = sum(sample["data"])
                data_card.exp[channel["name"]].update({sample["name"]: data})

                mods = [
                    spec["channels"][idxc]["samples"][idxs]["modifiers"][i]["type"]
                    for i, mod in enumerate(
                        spec["channels"][idxc]["samples"][idxs]["modifiers"]
                    )
                ]
                if "histosys" in mods:
                    # if shape uncertainty histogram is available
                    data_card.shapeMap[channel["name"]].update(
                        {
                            spec["channels"][idxc]["samples"][idxs]["name"]: [
                                shapefile,
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
                    data_card.shapeMap[channel["name"]].update(
                        {
                            spec["channels"][idxc]["samples"][idxs]["name"]: [
                                shapefile,
                                channel["name"]
                                + "/"
                                + spec["channels"][idxc]["samples"][idxs]["name"],
                            ]
                        }
                    )
                if "staterror" in mods:
                    # if staterror is present, provide variances for histograms
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
                    data_card.binParFlags.update({channel["name"]: True})
                elif "shapesys" in mods:  # for shapesys, also provide variances
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
                    data_card.binParFlags.update({channel["name"]: True})
                else:  # no staterror or shapesys->0 variance
                    h_data = hist.Hist.new.Regular(
                        channel_bins[channel["name"]], 0, channel_bins[channel["name"]]
                    ).Weight()
                    h_data[...] = np.stack(
                        [
                            spec["channels"][idxc]["samples"][idxs]["data"],
                            [0 for _ in range(channel_bins[channel["name"]])],
                        ],
                        axis=-1,
                    )
                file[
                    channel["name"]
                    + "/"
                    + spec["channels"][idxc]["samples"][idxs]["name"]
                ] = h_data


def addMods(file, spec, data_card, channel_bins, systs):
    """
    Add systematics to data_card
    """
    for im, modifier in enumerate(systs):
        if "normsys" in modifier[1]:  ##normsys
            # write normsys as 'shape?' so that Combine doesn't try to combine normsys and histosys mods of the same name
            data_card.systs.append((modifier[0], False, "shape?", [], {}))
            for channel in spec["channels"]:
                data_card.systs[im][4].update({channel["name"]: {}})
                for sample in channel["samples"]:
                    for i in data_card.systs:
                        i[4][channel["name"]].update({sample["name"]: 0.0})

        if "lumi" in modifier[1]:  ##lumi
            # Write lumi as lnN since they act the same way on the model
            data_card.systs.append((modifier[0], False, "lnN", [], {}))
            for channel in spec["channels"]:
                data_card.systs[im][4].update({channel["name"]: {}})
                for sample in channel["samples"]:
                    for i in data_card.systs:
                        i[4][channel["name"]].update({sample["name"]: 0.0})

        if "histosys" in modifier[1]:  ##histosys
            data_card.systs.append((modifier[0], False, "shape", [], {}))
            for channel in spec["channels"]:
                data_card.systs[im][4].update({channel["name"]: {}})
                for sample in channel["samples"]:
                    for i in data_card.systs:
                        i[4][channel["name"]].update({sample["name"]: 0.0})

    for idxc, channel in enumerate(spec["channels"]):
        for idxs, sample in enumerate(channel["samples"]):
            mods = sample["modifiers"]
            names = [mod["name"] for mod in mods]
            for syst in data_card.systs:
                name = syst[0]
                if name in names:
                    syst_type = syst[2]
                    # if systematic name is a modifier for this sample
                    if "shape?" in syst_type:  ##normsys
                        for mod in mods:
                            if mod["type"] == "normsys" and mod["name"] == name:
                                if mod["data"]["lo"] == 0:
                                    # asymmetric lnN
                                    syst[4][channel["name"]].update(
                                        {
                                            sample["name"]: str(
                                                mod["data"]["lo"] + 1e-9
                                            )
                                            + "/"
                                            + str(mod["data"]["hi"])
                                        }
                                    )
                                elif mod["data"]["hi"] == 0:
                                    # asymmetric lnN
                                    syst[4][channel["name"]].update(
                                        {
                                            sample["name"]: str(mod["data"]["lo"])
                                            + "/"
                                            + str(mod["data"]["hi"] + 1e-9)
                                        }
                                    )
                                else:
                                    # asymmetric lnN
                                    syst[4][channel["name"]].update(
                                        {
                                            sample["name"]: str(mod["data"]["lo"])
                                            + "/"
                                            + str(mod["data"]["hi"])
                                        }
                                    )
                    if "lnN" in syst_type:  ##lumi only
                        for mod in mods:
                            if mod["type"] == "lumi" and mod["name"] == name:
                                for measurement in spec["measurements"]:
                                    for param in measurement["config"]["parameters"]:
                                        if mod["name"] == param["name"]:
                                            # asymmetric lnN
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
                                                }
                                            )

                    if "shape" in syst_type:  ##histosys
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
                                            for _ in range(
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
                                            for _ in range(
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


def addSignal(spec, data_card, channels, modifiers):
    """
    Determine which samples are signal
    """
    measurements = [
        measurement["config"]["poi"] for measurement in spec["measurements"]
    ]
    signal_mods = [modifier[0] for modifier in modifiers if modifier[0] in measurements]

    for idxc, _ in enumerate(channels):
        for idxs, sample in enumerate(spec["channels"][idxc]["samples"]):
            for mod in spec["channels"][idxc]["samples"][idxs]["modifiers"]:
                for sig in signal_mods:
                    if sig == mod["name"]:
                        data_card.isSignal.update({sample["name"]: True})
                        data_card.signals.append(sample["name"])


def addRateParams(spec, data_card, channels, modifiers):
    """
    Add normfactor mods as rateParams (excluding signal strength).
    """
    measurements = [
        measurement["config"]["poi"] for measurement in spec["measurements"]
    ]
    signal_mods = [modifier[0] for modifier in modifiers if modifier[0] in measurements]

    for idxc, channel in enumerate(channels):
        for idxs, sample in enumerate((spec["channels"][idxc]["samples"])):
            is_signal = any(mod["name"] in signal_mods for mod in sample["modifiers"])
            if not is_signal:
                for mod in spec["channels"][idxc]["samples"][idxs]["modifiers"]:
                    # normfactor or shapefactor
                    if "normfactor" in mod["type"] or "shapefactor" in mod["type"]:
                        for measurement in spec["measurements"]:
                            for param in measurement["config"]["parameters"]:
                                data_card.rateParams.update(
                                    {f"{channel}AND" + sample["name"]: []}
                                )
                                if mod["name"] == param["name"]:
                                    data_card.rateParams[
                                        f"{channel}AND" + sample["name"]
                                    ].append([[mod["name"], 1, 0, param["bounds"]], ""])
                                else:
                                    data_card.rateParams[
                                        f"{channel}AND" + sample["name"]
                                    ].append([[mod["name"], 1, 0], ""])


def write_data_card(spec, data_card, channels, path):
    """
    Manually write each line of the datacard from data_card object
    """
    with open(path, "w") as f:
        f.write(f"imax {str(size(data_card.bins))}" + "\n")
        f.write(
            "jmax "
            + str(size(data_card.processes) - size(data_card.isSignal.keys()))
            + "\n"
        )
        f.write(f"kmax {str(size(data_card.systs, 0))}" + "\n")

        if data_card.hasShapes:
            for channel in data_card.shapeMap.keys():
                for sample in data_card.shapeMap[channel].keys():
                    f.write(
                        f"shapes {sample}  {channel}  {data_card.shapeMap[channel][sample][0]}  {data_card.shapeMap[channel][sample][1]}"
                    )
                    if size(data_card.shapeMap[channel][sample]) > 2:
                        f.write(f"  {data_card.shapeMap[channel][sample][2]}" + "\n")
                    else:
                        f.write("\n")

        f.write("\n---------------------------------\n")
        f.write("bin ")
        for bin in data_card.obs.keys():
            f.write(f"{bin} ")
        f.write("\n")
        f.write("observation ")
        for channel in data_card.obs.keys():
            f.write(f"{str(data_card.obs[channel])} ")
        f.write("\n---------------------------------\n")
        f.write("bin     ")
        for channel in data_card.obs.keys():
            for sample in data_card.exp[channel].keys():
                f.write(f"{channel}    ")
        f.write("\n")
        f.write("process     ")
        for channel in data_card.bins:
            for sample in data_card.exp[channel].keys():
                f.write(f"{sample}    ")
        f.write("\n")
        f.write("process     ")
        for channel in data_card.bins:
            for sample in data_card.exp[channel].keys():
                if sample in data_card.signals:
                    f.write(f"{str(-1 * data_card.processes.index(sample))}     ")
                else:
                    f.write(f"{str(data_card.processes.index(sample) + 1)}     ")
        f.write("\n")
        f.write("rate     ")
        for channel in data_card.bins:
            for sample in data_card.exp[channel].keys():

                f.write(f"{str(data_card.exp[channel][sample])}     ")
        f.write("\n---------------------------------\n")
        for syst in data_card.systs:
            f.write(f"{syst[0]}  {syst[2]}  ")
            for bin in syst[4].keys():
                for sample in data_card.exp[bin].keys():
                    if syst[4][bin][sample] != 0:
                        f.write(f"{str(syst[4][bin][sample])}  ")
                    else:
                        f.write("-       ")

            f.write("\n")
        f.write("\n---------------------------------\n")
        for cAp in data_card.rateParams.keys():
            _dir = cAp.split("AND")
            for i in range(size(data_card.rateParams[cAp], 0)):
                if size(data_card.rateParams[cAp][i][0]) > 3:
                    f.write(
                        f"{str(data_card.rateParams[cAp][i][0][0])} rateParam {_dir[0]} {_dir[1]} {str(data_card.rateParams[cAp][i][0][1])} {data_card.rateParams[cAp][i][0][3]}"
                    )
                else:
                    f.write(
                        f"{str(data_card.rateParams[cAp][i][0][0])} rateParam {_dir[0]} {_dir[1]} {str(data_card.rateParams[cAp][i][0][1])}"
                    )
                f.write("\n")
        f.write("\n---------------------------------\n")
        for idxc, channel in enumerate(channels):
            if (
                channel in data_card.binParFlags.keys()
                and data_card.binParFlags[channel] == True
            ):
                # double check to be safe
                shapesys = False
                staterror = False
                for sample in spec["channels"][idxc]["samples"]:
                    mod_types = [mod["type"] for mod in sample["modifiers"]]
                    if "shapesys" in mod_types:
                        shapesys = True
                    elif "staterror" in mod_types:
                        staterror = True

                if shapesys:
                    f.write(f"{channel} autoMCStats 100000 0 2" + "\n")
                if staterror:
                    f.write(f"{channel} autoMCStats 0 0 2" + "\n")


def pyhf_convert_to_datacard(workspace, outdatacard, shapefile):
    with open(workspace) as serialized:
        spec = json.load(serialized)
    workspace = pyhf.Workspace(spec)
    model = workspace.model()

    channels = model.config.channels
    channel_bins = model.config.channel_nbins
    samples = model.config.samples
    modifiers = model.config.modifiers
    lumi = model.config.parameters

    # generate list of mods without normfactors, shapesys, staterror
    systs = [
        mod
        for mod in modifiers
        if "normfactor" not in mod and "shapesys" not in mod and "staterror" not in mod
    ]

    data_card = Datacard()

    file = uproot.recreate(shapefile)
    addChannels(file, spec, data_card, channel_bins)
    addSamples(file, spec, data_card, channel_bins, samples, shapefile)
    addMods(file, spec, data_card, channel_bins, systs)
    addSignal(spec, data_card, channels, modifiers)
    addRateParams(spec, data_card, channels, modifiers)
    file.close()
    write_data_card(spec, data_card, channels, outdatacard)


def main():
    # Add command line args
    parser = OptionParser()
    parser.add_option(
        "-O",
        "--out-datacard",
        dest="outdatacard",
        default="converted_datacard.txt",
        help="desired name of datacard file",
    )
    parser.add_option(
        "-s",
        "--shape-file",
        dest="shapefile",
        default="shapes.root",
        help="desired name of shapes file",
    )
    options, args = parser.parse_args()

    pyhf_convert_to_datacard(
        workspace=args[0], outdatacard=options.outdatacard, shapefile=options.shapefile
    )


if __name__ == "__main__":
    main()
