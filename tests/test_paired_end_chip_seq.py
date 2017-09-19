import pytest

from utils.helpers import fetch_main_targets
from .helpers import run_dag


targets = fetch_main_targets()

@pytest.mark.dryrun
@pytest.mark.parametrize("target", targets)
def test_dna_repair_dag(target):

    exit_status = run_dag(target, "tests/test_data/paired_end/config.yaml")

    assert exit_status == 0


@pytest.mark.dryrun
@pytest.mark.bam
@pytest.mark.parametrize("target", targets)
def test_dna_repair_dag_bam(target):

    exit_status = run_dag(target, "tests/test_data/paired_end/config.yaml",
                          extras="--config filetype=bam sample_sheet=tests/test_data/paired_end/sample_sheet_bam.txt")

    assert exit_status == 0
