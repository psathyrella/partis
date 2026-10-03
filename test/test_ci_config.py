import os
import subprocess

import yaml

from _helpers import REPO_DIR

WORKFLOW = os.path.join(REPO_DIR, '.github', 'workflows', 'suite.yml')


def load_workflow():
    with open(WORKFLOW) as wfile:
        return yaml.safe_load(wfile)


def run_cmds(wflow):
    return [step.get('run', '') for job in wflow['jobs'].values() for step in job['steps']]


def test_triggers_only_on_prs_into_feature_branch():
    wflow = load_workflow()
    triggers = wflow.get('on', wflow.get(True))  # yaml 1.1 reads a bare 'on' key as True
    assert set(triggers) == {'pull_request'}
    assert triggers['pull_request']['branches'] == ['disjoint-grouping']


def test_installs_partis_and_runs_whole_suite():
    cmds = run_cmds(load_workflow())
    assert any('pip install -e .' in cmd for cmd in cmds)  # builds the c++ and installs xxhash from setup.py
    assert any('import xxhash' in cmd for cmd in cmds)
    assert any(cmd.strip() == 'make test' for cmd in cmds)


def test_make_fails_when_a_pipeline_script_fails(tmp_path):
    def make_pipelines(scripts):
        return subprocess.run(['make', '-C', REPO_DIR, 'pipeline-tests', 'PIPELINE_TESTS=%s' % ' '.join(scripts)], capture_output=True, text=True)
    passing, failing = tmp_path / 'test_pass_pipeline.sh', tmp_path / 'test_fail_pipeline.sh'
    passing.write_text('exit 0\n')
    failing.write_text('exit 1\n')
    assert make_pipelines([str(passing)]).returncode == 0
    assert make_pipelines([str(failing), str(passing)]).returncode != 0  # failure not masked by a later pass
