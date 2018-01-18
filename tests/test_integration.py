
import pytest

from utils.helpers import fetch_main_targets, multi_group_targets
from .helpers import run_dag


targets = fetch_main_targets()


@pytest.mark.dryrun
@pytest.mark.parametrize("target", targets + multi_group_targets)
def test_integration_dag(target):

    exit_status = run_dag(target, "tests/test_data/integration/config.yaml", "tests/test_data/integration/sample_sheet.txt", dryrun=False)

    assert exit_status == 0


# @pytest.mark.dryrun
# @pytest.mark.bam
# @pytest.mark.parametrize("target", targets + multi_group_targets)
# def test_integration_dag_bam(target):

#     exit_status = run_dag(target, "tests/test_data/integration/config.yaml",
#                           "tests/test_data/integration/sample_sheet_bam.txt",
#                           extras="--config filetype=bam")

#     assert exit_status == 0



# region_targets = []
# for target in targets + multi_group_targets:
#     if "heatmap" in target or "profileplot" in target:
#         region_targets.append(target)

# @pytest.mark.dryrun
# @pytest.mark.bam
# @pytest.mark.parametrize("target", region_targets)
# def test_integration_dag_bam(target):

#     exit_status = run_dag(target, "tests/test_data/integration/config.yaml",
#                           "tests/test_data/integration/sample_sheet.txt",
#                           extras="", dryrun=False)

#     assert exit_status == 0
