"""Checks for cohort identity, diagnostic pairing and false convergence claims."""
import copy
from pathlib import Path
import tempfile
import unittest

import collect_resolution_4097_results as collect
import report_resolution_4097 as report

ROOT = Path(__file__).resolve().parent.parents[4]/'oxygen/results/eFAST/neighborhood10pct/pilot/N4097'


class N4097ReportChecks(unittest.TestCase):
    def test_manifest_rejects_wrong_endpoint_or_focal_parameter(self):
        rows = collect.read(ROOT/'slurm/task_manifest.tsv')
        collect.validate_tasks(rows)
        bad = copy.deepcopy(rows);bad[0]['fit_seed']='seed101'
        with self.assertRaises(AssertionError): collect.validate_tasks(bad)
        bad = copy.deepcopy(rows);bad[1]['focal_parameter']='mu_hp'
        with self.assertRaises(AssertionError): collect.validate_tasks(bad)
        bad = copy.deepcopy(rows);bad[-1]['phase_seed']=rows[0]['phase_seed']
        with self.assertRaises(AssertionError): collect.validate_tasks(bad)

    def test_primary_validation_rejects_the_old_primary_pair_label(self):
        _,values,_,oxygen = report.load(ROOT)
        rows = report.read(ROOT/'resolution_diagnostics.tsv')
        changed = next(r for r in rows if int(r['N_low'])==2049)
        self.assertEqual(changed['primary_comparison'],'True')
        changed['primary_comparison']='False'
        with tempfile.TemporaryDirectory() as directory:
            temp = Path(directory)
            for p in ROOT.iterdir(): (temp/p.name).symlink_to(p,target_is_directory=p.is_dir())
            (temp/'resolution_diagnostics.tsv').unlink()
            report.write(temp/'resolution_diagnostics.tsv',rows)
            with self.assertRaisesRegex(AssertionError,'primary_comparison'):
                report.validate(temp,values,oxygen)

    def test_joint_pass_uses_the_new_resolution_and_preserves_local_failure(self):
        repeats = report.read(ROOT/'repeat_diagnostics.tsv')
        deltas = report.read(ROOT/'resolution_diagnostics.tsv')
        for row in repeats:
            if row['fit_seed']=='seed25' and row['parameter']=='buffer_beta' and row['index']=='ST' and row['output']=='dominant_mean_ploidy' and row['oxygen_window']=='overall_0_to_5pct':
                if int(row['N'])==2049: row.update(p90_pass='False',diagnostic_status='not_passed')
                if int(row['N'])==4097: row.update(p90_pass='True',diagnostic_status='localized_review',phase_range_p90='.09',phase_range_max='.50')
        for row in deltas:
            if row['fit_seed']=='seed25' and row['parameter']=='buffer_beta' and row['index']=='ST' and row['output']=='dominant_mean_ploidy' and row['oxygen_window']=='overall_0_to_5pct' and int(row['N_low'])==2049:
                row.update(p90_pass='True',diagnostic_status='localized_review',absolute_delta_p90='.01',absolute_delta_max='.20')
        rows = report.detail_rows(repeats,deltas)
        item = next(r for r in rows if r['fit_seed']=='seed25' and r['parameter']=='buffer_beta' and r['oxygen_window']=='overall_0_to_5pct')
        self.assertTrue(item['N4097_both_p90_pass'])
        self.assertFalse(item['N4097_both_max_pass'])
        self.assertEqual(len(rows),6*14*3)


if __name__ == '__main__': unittest.main()
