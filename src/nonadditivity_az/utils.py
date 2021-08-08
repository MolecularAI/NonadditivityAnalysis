import multiprocessing as mp
from pathlib import Path
from typing import Tuple

import numpy as np
import pandas as pd
import pystow
from rdkit import Chem
from rdkit.Chem import MolStandardize
from rdkit.Chem.MolStandardize import rdMolStandardize

import chembl_downloader
from chembl_downloader.queries import get_assay_sql

MODULE = pystow.module("nonadditivity")


def get_processed_assay_df(assay_chembl_id: str) -> Tuple[pd.DataFrame, Path]:
    submodule = MODULE.submodule(assay_chembl_id)
    assay_raw_path = submodule.join(name="raw.tsv")
    assay_processed_path = submodule.join(name="processed.tsv")

    if assay_processed_path.is_file():
        return pd.read_csv(assay_processed_path, sep="\t"), assay_processed_path

    if assay_raw_path.is_file():
        data = pd.read_csv(assay_raw_path, sep="\t")
    else:
        data = chembl_downloader.query(get_assay_sql(assay_chembl_id))
        data.columns = [
            "Smiles",
            "Molecule ChEMBL ID",
            "Standard Type",
            "Standard Relation",
            "Standard Value",
            "Standard Units",
        ]
        data.to_csv(assay_raw_path, sep="\t", index=False)

    df = process(data)
    df.to_csv(assay_processed_path, sep="\t", index=False)
    return df, assay_processed_path


def discard_nan_smiles(df):
    """Remove NaNs from SMILES column"""
    df = df.dropna(subset=["SMILES"])
    return df


def discard_uncertain_values(df):
    """Deleting uncertain values"""
    df = df[df["Standard Relation"] != "'>'"]
    df = df[df["Standard Relation"] != "'<'"]
    df["VALUE"] = df["VALUE"].astype(float)
    df = df[df["VALUE"] > 0]  # in case one needs to delete negative values
    df = df.drop(columns=["Standard Relation"])
    return df


unit_conversion = {
    "M": 1,
    "mM": 1000,
    "uM": 1000000,
    "nM": 1000000000,
    "pM": 1000000000000,
    "fM": 1000000000000000,
}


def log_converstion(x, UNIT):
    if UNIT not in unit_conversion:
        return x
    x = -1 * np.log10(x / unit_conversion[UNIT])
    return x


def create_conversion_column(df):
    """Converting values to logged ones"""
    arr = [log_converstion(x["VALUE"], x["UNIT"]) for idx, x in df.iterrows()]

    df["NEW_VALUE"] = arr
    df = df[df["NEW_VALUE"] > 0]
    df = df.drop(columns=["VALUE"])

    # deleting the values that are more than 10 mM and less than 1 fM

    df = df[df["NEW_VALUE"] > 2]
    df = df[df["NEW_VALUE"] < 11]

    return df


def calculate_average(df):
    """Calculating the avarage of the activity and calculating the median"""
    df["Median_Value"] = df.groupby(["COMPOUND_NAME"])["NEW_VALUE"].transform("median")
    df["max_value"] = df.groupby(["COMPOUND_NAME"])["NEW_VALUE"].transform("max")
    df["min_value"] = df.groupby(["COMPOUND_NAME"])["NEW_VALUE"].transform("min")
    df["difference"] = df.max_value - df.min_value

    df = df.drop_duplicates(subset=["COMPOUND_NAME"], keep="first")

    return df


def discard_ambiguous_compound_measurements(df, max_thrs=2.5):
    """Delete compounds that have been measured several times in one test and differ more than 2.5 log units"""
    df = df[df.difference < max_thrs]

    df = df[["SMILES", "COMPOUND_NAME", "ENDPOINT", "Median_Value", "MEASUREMENT"]]
    df = df.rename(columns=({"Median_Value": "VALUE"}))

    return df


def standardize_rdkit(row, col):
    """Standardize molecules using RDkit"""
    smi = row[col]

    try:
        mol = Chem.MolFromSmiles(smi)  # sanitization is done by default
        fmol = rdMolStandardize.FragmentParent(mol)  # returns largest fragment
        cmol = rdMolStandardize.ChargeParent(fmol)  # uncharges the largest fragment
        smi = Chem.MolToSmiles(cmol)
        ssmi = MolStandardize.canonicalize_tautomer_smiles(
            smi
        )  # returns the canonicalized tautomer
        tsmi = MolStandardize.rdMolStandardize.StandardizeSmiles(ssmi)  # standardize
    except:
        tsmi = "none"

    return tsmi


def generateStandarizedSmiles(smilesfile, smiles_column):
    # set number of cores for parallelization
    pool = mp.Pool(8)
    stsmi_list = pool.starmap(
        standardize_rdkit, [(smi, smiles_column) for idx, smi in smilesfile.iterrows()]
    )
    pool.close()

    smilesfile[smilesfile.columns[smiles_column]] = stsmi_list

    return smilesfile


def merge_duplicate_smiles(df):
    """Discard duplicate SMILES

    - Keep the one with the highest value, i.e. most active one
    """
    df = df.sort_values("VALUE").drop_duplicates(subset=["SMILES"], keep="last")
    return df


def discarding_heavy_mols(smi, min_size=0, max_size=70):
    """Remove molecules with > 70 heavy atoms"""
    try:
        mol = Chem.MolFromSmiles(smi, sanitize=False)
        if min_size <= mol.GetNumHeavyAtoms() <= max_size:
            return False
        else:
            return True
    except Exception:
        return True


def removeHeavyMols(df, smiles_column):
    idx = []
    discard = []
    for index, row in df.iterrows():
        # for smi in df.iloc[:,smiles_column]:
        if discarding_heavy_mols(row[smiles_column]):
            idx.append(index)
            discard.append(row.values.tolist())

    df.drop(idx, inplace=True)
    return df


def process(data: pd.DataFrame):
    df = data.rename(
        columns=(
            {
                "Molecule ChEMBL ID": "COMPOUND_NAME",
                "Smiles": "SMILES",
                "Standard Value": "VALUE",
                "Standard Units": "UNIT",
                "Standard Type": "ENDPOINT",
            }
        )
    )
    df = df[
        ["SMILES", "COMPOUND_NAME", "ENDPOINT", "Standard Relation", "VALUE", "UNIT"]
    ]

    print("#cmpds: ", len(df["COMPOUND_NAME"]))
    print("#unique cmpds: ", len(df["COMPOUND_NAME"].value_counts()))

    # Counting how many times compounds were measured in tests
    df["MEASUREMENT"] = df.groupby(["COMPOUND_NAME"])["COMPOUND_NAME"].transform(
        "count"
    )

    # Discard cmpds without SMILES
    df = discard_nan_smiles(df)
    print("#cpds with SMILES: ", len(df.iloc[:, 0]))

    # Discard ambiguous data
    df = discard_uncertain_values(df)
    print("#cpds with values: ", len(df.iloc[:, 0]))

    # Convert IC50 to pIC50
    df = create_conversion_column(df)

    # Calculate average values and discard cpds with > 2.5 log unit measurement differences
    df = calculate_average(df)
    df = discard_ambiguous_compound_measurements(df)
    print("#cpds after merging multi measurements: ", len(df.iloc[:, 0]))

    # standardize SMILES, merge duplicates and retain higher active one
    smiles_column = 0
    df = generateStandarizedSmiles(df, smiles_column)
    df = df[df["SMILES"] != "none"]
    df = merge_duplicate_smiles(df)
    print("#cpds after merging duplicate SMILES: ", len(df.iloc[:, 0]))

    # Remove cpds with > 70 HA
    smiles_column = 0
    df = removeHeavyMols(df, smiles_column)
    print("#cpds < 70 HA: ", len(df.iloc[:, 0]))

    # Rename columns for subsequent NAA
    df = df.rename(columns=({"COMPOUND_NAME": "ID"}))
    df = df[["ID", "SMILES", "VALUE", "MEASUREMENT"]]
    return df
