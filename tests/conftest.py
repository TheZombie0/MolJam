import os
import pytest
import pandas as pd

FIXTURES_DIR = os.path.join(os.path.dirname(__file__), "fixtures")


@pytest.fixture
def tiny_db_path():
    return os.path.join(FIXTURES_DIR, "tiny_db.csv")


@pytest.fixture
def tiny_db_df(tiny_db_path):
    return pd.read_csv(tiny_db_path)


@pytest.fixture
def scorer(tiny_db_df):
    from moljam.scoring.scorer import MoleculeDBScorer

    return MoleculeDBScorer(
        tiny_db_df,
        smiles_col="smiles",
        activity_cols=["activity"],
        class_cols=["label"],
        experimental_method_cols=["exp_method"],
        id_col="mol_id",
        name_col="mol_name",
        time_col="time_point",
        use_parallel=False,
        include_experimental_info=False,
    )


@pytest.fixture
def scorer_with_experimental(tiny_db_df):
    from moljam.scoring.scorer import MoleculeDBScorer

    return MoleculeDBScorer(
        tiny_db_df,
        smiles_col="smiles",
        activity_cols=["activity"],
        class_cols=["label"],
        experimental_method_cols=["exp_method"],
        id_col="mol_id",
        name_col="mol_name",
        time_col="time_point",
        use_parallel=False,
        include_experimental_info=True,
    )
