import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from PIL import Image, ImageDraw
from rdkit import Chem, Geometry
from rdkit.Chem import AllChem, Draw, rdFMCS
from scipy.stats import normaltest


def NA_distribution(x, significant_thrs, alpha: float = 0.05):
    sns.set_style("ticks")
    fig, ax = plt.subplots(1, 1, figsize=(6, 3.5))
    sns.distplot(
        x,
        hist=True,
        kde=True,
        bins=int(100 / 5),
        color="crimson",
        kde_kws={"shade": True, "linewidth": 2},
        ax=ax,
        label="Real",
    )
    n_of_cycles = len(x)
    normal_dist = np.random.normal(loc=0, scale=significant_thrs, size=n_of_cycles)
    sns.distplot(
        normal_dist,
        hist=False,
        kde=True,
        color="grey",
        kde_kws={"shade": True, "linewidth": 2},
        ax=ax,
        label="Theoretical",
    )
    plt.legend()

    stat, p = normaltest(x)
    if p > alpha:
        plt.title(f"Nonadditivity is Gaussian ($p<{p:.2e}$)")
    else:
        plt.title(f"Nonadditivity is not Gaussian ($p<{p:.2e}$)")

    plt.xlabel("Nonadditivity")
    plt.ylabel("Density")

    return fig, ax


def plot_outliers(x, title):
    sns.set_style("ticks")
    fig, ax = plt.subplots(1, 1)
    sns.distplot(
        x,
        hist=True,
        kde=True,
        color="crimson",
        kde_kws={"shade": True, "linewidth": 2},
        ax=ax,
    )

    plt.title(title)
    plt.xlabel("pActivity Value")
    plt.ylabel("Density")


def draw_image(
    ids,
    smiles,
    tsmarts,
    pActs,
    Acts,
    qualifiers,
    nonadd,
    mcss_tot,
    image_file,
    target="",
):
    """
    Draw Nonadditivity Circle to Image file
    """

    cpds = [Chem.MolFromSmiles(i) for i in smiles]

    #########
    # Compute Coordinates of local MCSS, aligned with global MCSS

    mcss_loc = Chem.MolFromSmarts(
        rdFMCS.FindMCS(cpds, completeRingsOnly=True, timeout=60).smartsString
    )
    Chem.GetSymmSSSR(mcss_loc)

    if mcss_tot:
        mcss_tot_coords = [
            mcss_tot.GetConformer().GetAtomPosition(x)
            for x in range(mcss_tot.GetNumAtoms())
        ]
        coords2D_tot = [Geometry.Point2D(pt.x, pt.y) for pt in mcss_tot_coords]

        mcss_match = mcss_loc.GetSubstructMatch(mcss_tot)

        coordDict = {}
        for i, coord in enumerate(coords2D_tot):
            coordDict[mcss_match[i]] = coord

        AllChem.Compute2DCoords(mcss_loc, coordMap=coordDict)
    else:
        AllChem.Compute2DCoords(mcss_loc)

    #########
    # Align circle compounds to local MCSS

    matchVs = [x.GetSubstructMatch(mcss_loc) for x in cpds]

    # compute reference coordinates:
    mcss_loc_coords = [
        mcss_loc.GetConformer().GetAtomPosition(x)
        for x in range(mcss_loc.GetNumAtoms())
    ]
    coords2D_loc = [Geometry.Point2D(pt.x, pt.y) for pt in mcss_loc_coords]

    # generate coords for the other molecules using the common substructure
    for molIdx in range(4):
        mol = cpds[molIdx]
        coordDict = {}
        for i, coord in enumerate(coords2D_loc):
            coordDict[matchVs[molIdx][i]] = coord
        AllChem.Compute2DCoords(mol, coordMap=coordDict)

    ##########
    # Assemble Image

    qualifiers_inv = ["" for i in range(4)]
    for i in range(4):
        if qualifiers[i] == ">":
            qualifiers_inv[i] = "<"
        elif qualifiers[i] == "<":
            qualifiers_inv[i] = ">"
        else:
            continue

    new_im = Image.new("RGB", size=(650, 670), color=(255, 255, 255, 0))
    if Acts[0] != "":
        new_im.paste(
            Draw.MolToImage(
                cpds[0],
                size=(300, 300),
                legend=ids[0]
                + "        "
                + qualifiers_inv[0]
                + Acts[0]
                + " ("
                + qualifiers[0]
                + pActs[0]
                + ")",
            ),
            (0, 0),
        )
        new_im.paste(
            Draw.MolToImage(
                cpds[1],
                size=(300, 300),
                legend=ids[1]
                + "        "
                + qualifiers_inv[1]
                + Acts[1]
                + " ("
                + qualifiers[1]
                + pActs[1]
                + ")",
            ),
            (350, 0),
        )
        new_im.paste(
            Draw.MolToImage(
                cpds[2],
                size=(300, 300),
                legend=ids[2]
                + "        "
                + qualifiers_inv[2]
                + Acts[2]
                + " ("
                + qualifiers[2]
                + pActs[2]
                + ")",
            ),
            (350, 350),
        )
        new_im.paste(
            Draw.MolToImage(
                cpds[3],
                size=(300, 300),
                legend=ids[3]
                + "        "
                + qualifiers_inv[3]
                + Acts[3]
                + " ("
                + qualifiers[3]
                + pActs[3]
                + ")",
            ),
            (0, 350),
        )

        draw = ImageDraw.Draw(new_im)
        # font = ImageFont.truetype(font_path, 14)
        draw.text(
            (260, 330), "Nonadditivity: " + nonadd, fill=(0, 0, 0, 0)
        )  # , font=font)

        # font = ImageFont.truetype(font_path, 10)
        if target != "":
            draw.text(
                (10, 650),
                "[uM]  (-log10[M])  Activity in Assay: " + target,
                fill=(0, 0, 0, 0),
            )  # , font=font)
    else:
        new_im.paste(
            Draw.MolToImage(
                cpds[0],
                size=(300, 300),
                legend=ids[0] + "        " + qualifiers[0] + pActs[0],
            ),
            (0, 0),
        )
        new_im.paste(
            Draw.MolToImage(
                cpds[1],
                size=(300, 300),
                legend=ids[1] + "        " + qualifiers[1] + pActs[1],
            ),
            (350, 0),
        )
        new_im.paste(
            Draw.MolToImage(
                cpds[2],
                size=(300, 300),
                legend=ids[2] + "        " + qualifiers[2] + pActs[2],
            ),
            (350, 350),
        )
        new_im.paste(
            Draw.MolToImage(
                cpds[3],
                size=(300, 300),
                legend=ids[3] + "        " + qualifiers[3] + pActs[3],
            ),
            (0, 350),
        )

        draw = ImageDraw.Draw(new_im)
        # font = ImageFont.truetype(font_path, 14)
        draw.text(
            (260, 330), "Nonadditivity: " + nonadd, fill=(0, 0, 0, 0)
        )  # , font=font)

        # font = ImageFont.truetype(font_path, 10)
        if target != "":
            draw.text(
                (10, 650), "Activity in Assay: " + target, fill=(0, 0, 0, 0)
            )  # , font=font)

    # Draw Arrows
    draw.line((300, 150, 350, 150), fill=0, width=2)
    draw.line((340, 145, 350, 150), fill=0, width=2)
    draw.line((340, 155, 350, 150), fill=0, width=2)

    draw.line((300, 500, 350, 500), fill=0, width=2)
    draw.line((340, 495, 350, 500), fill=0, width=2)
    draw.line((340, 505, 350, 500), fill=0, width=2)

    draw.line((150, 300, 150, 350), fill=0, width=2)
    draw.line((145, 340, 150, 350), fill=0, width=2)
    draw.line((155, 340, 150, 350), fill=0, width=2)

    draw.line((500, 300, 500, 350), fill=0, width=2)
    draw.line((495, 340, 500, 350), fill=0, width=2)
    draw.line((505, 340, 500, 350), fill=0, width=2)

    # Add Reaction Parts
    b = Chem.MolFromSmiles(tsmarts[0][: tsmarts[0].index(">")])
    new_im.paste(Draw.MolToImage(b, size=(50, 50)), (300, 90))

    b = Chem.MolFromSmiles(tsmarts[0][tsmarts[0].index(">") + 2 :])
    new_im.paste(Draw.MolToImage(b, size=(50, 50)), (300, 160))

    b = Chem.MolFromSmiles(tsmarts[0][: tsmarts[0].index(">")])
    new_im.paste(Draw.MolToImage(b, size=(50, 50)), (300, 440))

    b = Chem.MolFromSmiles(tsmarts[0][tsmarts[0].index(">") + 2 :])
    new_im.paste(Draw.MolToImage(b, size=(50, 50)), (300, 510))

    b = Chem.MolFromSmiles(tsmarts[1][: tsmarts[1].index(">")])
    new_im.paste(Draw.MolToImage(b, size=(50, 50)), (80, 300))

    b = Chem.MolFromSmiles(tsmarts[1][tsmarts[1].index(">") + 2 :])
    new_im.paste(Draw.MolToImage(b, size=(50, 50)), (170, 300))

    b = Chem.MolFromSmiles(tsmarts[1][: tsmarts[1].index(">")])
    new_im.paste(Draw.MolToImage(b, size=(50, 50)), (430, 300))

    b = Chem.MolFromSmiles(tsmarts[1][tsmarts[1].index(">") + 2 :])
    new_im.paste(Draw.MolToImage(b, size=(50, 50)), (520, 300))

    new_im.save(image_file)
