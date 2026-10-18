"""Guard scope, comparison labels, reference stability and localized failures."""
import copy
from pathlib import Path
import tempfile
import unittest

import numpy as np
import collect_resolution_nm_results as collect
import report_resolution_nm as report

ROOT=Path(__file__).resolve().parent.parents[4]/'oxygen/results/eFAST/neighborhood10pct/pilot'/collect.PILOT


class NMReportChecks(unittest.TestCase):
    def test_manifest_rejects_wrong_M_seed_or_phase(self):
        for config in collect.CONFIGS:
            rows=collect.read(ROOT/config/'slurm/task_manifest.tsv')
            collect.validate_tasks(rows,config)
            for field,value in [('M','4' if config=='N4097_M8' else '8'),('fit_seed','seed311'),('phase_seed','123')]:
                wrong=copy.deepcopy(rows);wrong[0][field]=value
                with self.assertRaises(AssertionError):collect.validate_tasks(wrong,config)

    def test_M_setting_comparison_rejects_N_resolution_label(self):
        config='N4097_M8'
        values,_,oxygen=report.load(ROOT,config)
        rows=report.read(ROOT/config/'configuration_comparison.tsv')
        self.assertEqual(rows[0]['comparison_kind'],'M_robustness')
        rows[0]['comparison_kind']='N_resolution'
        with tempfile.TemporaryDirectory() as directory:
            target=Path(directory)
            for path in (ROOT/config).iterdir():(target/path.name).symlink_to(path,target_is_directory=path.is_dir())
            (target/'configuration_comparison.tsv').unlink()
            report.write(target/'configuration_comparison.tsv',rows)
            with self.assertRaises(AssertionError):report.validate(target,values,oxygen,config)

    def test_test_setting_joint_pass_is_distinct_from_both_settings_repeat_pass(self):
        oxygen=np.arange(201)*.025
        for config,(n,m,seeds,kind) in collect.CONFIGS.items():
            values={}
            for seed in seeds:
                for level,harmonics in ((2049,4),(4097,4),(n,m)):
                    values[(seed,level,harmonics)]=np.full((5,2,14,201,2),.2)
                values[(seed,4097,4)][0,1,0,:,0]+=.2
            detail=report.detail_rows(config,values,oxygen)
            row=next(r for r in detail if r['fit_seed']==seeds[0] and r['parameter']==report.efast.ACTIVE[0]
                     and r['index']=='ST' and r['output']=='dominant_mean_ploidy' and r['oxygen_window']=='overall_0_to_5pct')
            self.assertTrue(row['test_repeat_and_difference_p90_pass'])
            self.assertFalse(row['both_configurations_R90_and_D90_pass'])
            self.assertFalse(row['R_reference_90_pass'])
            self.assertEqual((row['N_test'],row['M_test']),(n,m))
            if m==8:self.assertIsNone(row['delta_p90_improved_vs_previous_N_step'])
            values[(seeds[0],n,m)][0,1,0,77,0]+=1.1
            detail=report.detail_rows(config,values,oxygen)
            row=next(r for r in detail if r['fit_seed']==seeds[0] and r['parameter']==report.efast.ACTIVE[0]
                     and r['index']=='ST' and r['output']=='dominant_mean_ploidy' and r['oxygen_window']=='overall_0_to_5pct')
            self.assertTrue(row['test_repeat_and_difference_p90_pass'])
            self.assertFalse(row['test_repeat_and_difference_max_pass'])
            self.assertEqual(row['R_test_status'],'localized_review')
            self.assertEqual(row['D_status'],'localized_review')
            self.assertAlmostEqual(row['R_test_max_O2_pct'],1.925)


if __name__=='__main__':unittest.main()
