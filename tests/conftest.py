import os
import sys
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

TEST_DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
EXAMPLE_DIR = os.path.join(ROOT, "example")


def pytest_addoption(parser):
    parser.addoption(
        "--network",
        action="store_true",
        default=False,
        help="Run tests that require downloading data from Logan",
    )


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "network: requires network access (download from Logan); pass --network to enable",
    )


def pytest_runtest_setup(item):
    if "network" in item.keywords and not item.config.getoption("--network"):
        pytest.skip("network test — pass --network to run")


@pytest.fixture(scope="session")
def query_fa():
    return os.path.join(EXAMPLE_DIR, "query.fa")


@pytest.fixture(scope="session")
def accessions_txt():
    return os.path.join(EXAMPLE_DIR, "accessions.txt")


@pytest.fixture(scope="session")
def self_blast_txt():
    return os.path.join(TEST_DATA, "self_blast.txt")


@pytest.fixture(scope="session")
def expected_self_synth():
    with open(os.path.join(TEST_DATA, "expected_self_synth.txt")) as f:
        return f.read()
