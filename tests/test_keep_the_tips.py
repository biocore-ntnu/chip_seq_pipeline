import pytest

from utils.helpers import fetch_main_targets, multi_group_targets
from .helpers import run_dag


targets = fetch_main_targets()


@pytest.mark.dryrun
@pytest.mark.parametrize("target", targets + multi_group_targets)
def test_dna_repair_dag(target):

    exit_status = run_dag(target, "tests/test_data/keep_the_tips/config.yaml")

    assert exit_status == 0


@pytest.mark.dryrun
@pytest.mark.bam
@pytest.mark.parametrize("target", targets + multi_group_targets)
def test_keep_the_tips_dag_bam(target):

    exit_status = run_dag(target, "tests/test_data/keep_the_tips/config.yaml",
                          extras="--config bam=true sample_sheet=tests/test_data/keep_the_tips/sample_sheet.txt")

    assert exit_status == 0