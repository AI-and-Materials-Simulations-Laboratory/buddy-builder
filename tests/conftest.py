import os
from pathlib import Path

# Ask Buddy: never read a developer .env (it may hold a real key) and never call a model in tests.
os.environ["BUDDY_DOTENV"] = "0"
for _k in ("BUDDY_LLM_API_KEY", "BUDDY_LLM_PROVIDER", "OPENAI_API_KEY", "ANTHROPIC_API_KEY",
           "ATOMFORGE_LLM_API_KEY", "ATOMFORGE_LLM_PROVIDER"):
    os.environ.pop(_k, None)

import pytest

from nanobuilder.crystal import build_crystal
from nanobuilder.ligands import build_ligand
from nanobuilder.schema import GeometryConfig, LigandConfig, MaterialConfig
from nanobuilder.shapes import build_geometry

FIXTURES = Path(__file__).resolve().parents[1] / "examples" / "fixtures"
EXAMPLES_DIR = Path(__file__).resolve().parents[1] / "examples"


@pytest.fixture(scope="session")
def fixtures():
    return FIXTURES


@pytest.fixture(scope="session")
def au():
    return build_crystal(MaterialConfig())


@pytest.fixture(scope="session")
def au_sphere(au):
    return build_geometry(au, GeometryConfig())


@pytest.fixture(scope="session")
def au111_slab(au):
    return build_geometry(au, GeometryConfig(shape="slab", miller=(1, 1, 1)))


@pytest.fixture(scope="session")
def scco_thiolate():
    return build_ligand(LigandConfig(smiles="SCCO", anchor_atom_index=0,
                                     binding_template="Au_thiolate", surface_species=["Au"]))


def upload_config(fixtures_dir, name, file_format="cif"):
    path = fixtures_dir / name
    return MaterialConfig(source="uploaded_file", preset_id=None, file_format=file_format,
                          filename=path.name, uploaded_bytes=path.read_bytes())
