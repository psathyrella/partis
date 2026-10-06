import os
import subprocess

from _helpers import REPO_DIR, read_yaml

WORKFLOW = os.path.join(REPO_DIR, '.github', 'workflows', 'suite.yml')


def load_workflow():
    return read_yaml(WORKFLOW)


def run_cmds(wflow):
    return [step.get('run', '') for job in wflow['jobs'].values() for step in job['steps']]


def test_triggers_only_on_prs_into_feature_branch():
    wflow = load_workflow()
    triggers = wflow.get('on', wflow.get(True))  # yaml 1.1 reads a bare 'on' key as True
    assert set(triggers) == {'pull_request'}
    assert triggers['pull_request']['branches'] == ['disjoint-grouping']


def test_runs_on_open_and_on_label_not_per_push():
    wflow = load_workflow()
    triggers = wflow.get('on', wflow.get(True))
    assert set(triggers['pull_request']['types']) == {'opened', 'reopened', 'ready_for_review', 'labeled'}
    assert wflow['jobs']['suite']['if'] == "github.event.action != 'labeled' || github.event.label.name == 'run-suite'"
    assert wflow['concurrency']['cancel-in-progress'] is True
    assert 'github.event.label.name' in wflow['concurrency']['group']


def test_installs_partis_and_runs_whole_suite():
    cmds = run_cmds(load_workflow())
    assert any('pip install -e .' in cmd for cmd in cmds)  # builds the c++ and installs xxhash from setup.py
    assert any('import xxhash' in cmd for cmd in cmds)
    assert any(cmd.strip() == 'make test' for cmd in cmds)


def test_runs_on_both_backends():
    job = load_workflow()['jobs']['suite']
    assert job['timeout-minutes'] > 0
    assert job['strategy']['matrix']['backend'] == ['cpp', 'zig']
    assert job['env']['PARTIS_TEST_BACKEND'] == '${{ matrix.backend }}'
    assert any(step.get('if') == "matrix.backend == 'zig'" and 'bin/zig-build.sh' in step['run'] for step in job['steps'])


def test_make_fails_when_a_pipeline_script_fails(tmp_path):
    def make_pipelines(scripts):
        return subprocess.run(['make', '-C', REPO_DIR, 'pipeline-tests', 'PIPELINE_TESTS=%s' % ' '.join(scripts)], capture_output=True, text=True)
    passing, failing = tmp_path / 'test_pass_pipeline.sh', tmp_path / 'test_fail_pipeline.sh'
    passing.write_text('exit 0\n')
    failing.write_text('exit 1\n')
    assert make_pipelines([str(passing)]).returncode == 0
    assert make_pipelines([str(failing), str(passing)]).returncode != 0  # failure not masked by a later pass
