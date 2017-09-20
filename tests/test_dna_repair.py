import pytest

from utils.helpers import fetch_main_targets, multi_group_targets
from .helpers import run_dag


targets = fetch_main_targets()
# targets = ["peaks", "input_profileplots"]


@pytest.mark.dryrun
@pytest.mark.parametrize("target", targets + multi_group_targets)
def test_dna_repair_dag(target):

    exit_status = run_dag(target, "tests/test_data/dna_repair/config.yaml", "tests/test_data/dna_repair/sample_sheet.txt")

    assert exit_status == 0


@pytest.mark.dryrun
@pytest.mark.bam
@pytest.mark.parametrize("target", targets + multi_group_targets)
def test_dna_repair_dag_bam(target):

    exit_status = run_dag(target, "tests/test_data/dna_repair/config.yaml",
                          "tests/test_data/dna_repair/sample_sheet_bam.txt",
                          extras="--config filetype=bam")

    assert exit_status == 0
