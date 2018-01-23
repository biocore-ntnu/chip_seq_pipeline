
def pytest_addoption(parser):
    parser.addoption("--configs", type=str,
        help="Additonal config flags for snakemake")

    parser.addoption("--cpu", type=int,
        help="Number of cores to use in integration tests.")

    parser.addoption("--targets", type=str, nargs="*",
        help="Makefile targets to test")

    parser.addoption("--forceall", action="store_true",
        help="Run Snakefile from very first job")

def pytest_generate_tests(metafunc):
    if 'integration_cpu' in metafunc.fixturenames:
        cores = metafunc.config.getoption('cpu') or 1
        metafunc.parametrize("integration_cpu", [cores])

    if 'configs' in metafunc.fixturenames:
        configs = metafunc.config.getoption('configs') or ""
        metafunc.parametrize("configs", [configs])

    # if 'targets' in metafunc.fixturenames:
    #     target = metafunc.config.getoption('targets')
    #     metafunc.parametrize("cli_target", target)

    if 'forceall' in metafunc.fixturenames:
        forceall = metafunc.config.getoption('forceall') or False
        metafunc.parametrize("forceall", [forceall])
