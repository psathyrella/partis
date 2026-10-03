import os

REPO_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
TEST_DIR = os.path.join(REPO_DIR, 'test')

# tracked paired simulation and the parameters inferred on it
PAIRED_SIMU_DIR = os.path.join(TEST_DIR, 'paired', 'ref-results', 'test', 'simu')
PAIRED_PARAM_DIR = os.path.join(TEST_DIR, 'paired', 'ref-results', 'test', 'parameters', 'simu')
