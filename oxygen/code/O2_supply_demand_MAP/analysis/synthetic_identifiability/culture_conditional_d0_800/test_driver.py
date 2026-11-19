"""Controller integration fixtures use no Slurm, model, or synthetic generation."""
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch
import driver
import workflow as w

class DriverTests(unittest.TestCase):
    def setUp(self):
        self.tmp=TemporaryDirectory();self.addCleanup(self.tmp.cleanup);self.s=Path(self.tmp.name)
        for name in ('state','tables','checks','engineering','data','design','truth'):(self.s/name).mkdir()
        self.rows=[];self.submitted=[];self.finished=[]
        for m in ('prepare','replay','parallel','engineering_check','replay_final','local_api_check','positive_local','positive_check'):self.add(m)
        for p in ['tables/engineering_gate.tsv','tables/search_controls.tsv','truth/manifest.tsv','data/expected_observations.rds','design/blind_bundle.rds']+['engineering/'+m+'/payload_manifest.tsv' for m in ('prepare','replay','parallel','engineering_check','replay_final')]:w.text(self.s/p,'fixture only\n')
        w.typed(self.s/'tables/engineering_gate.tsv',dict(**{'pass':True},cost_fits=True))
        w.typed(self.s/'checks/positive_check/result.tsv',dict(prediction_recovery=True,local_directional_stability=True))
        self.patches=[patch.object(w,'reconcile',side_effect=lambda s:self.rows),patch.object(w,'submit',side_effect=lambda s,t:self.submitted.append(t)),patch.object(driver,'accepted',return_value=None),patch.object(driver,'finish',side_effect=lambda s,r,e:self.finished.append(e))]
        for p in self.patches:p.start();self.addCleanup(p.stop)
    def add(self,mode,state='COMPLETED',task=None):
        self.rows.append(dict(task=task or mode,mode=mode,state=state,exit_code='0:0' if state=='COMPLETED' else '1:0',actual_core_seconds='1',reserved_core_seconds='100',phase='primary',job_id=str(len(self.rows)+1)))
    def chain(self,origin,good=False):
        for step in ('screen','global','globalcheck','local','check'):self.add(step+'_'+origin)
        w.typed(self.s/'checks'/('check_'+origin)/'result.tsv',dict(prediction_recovery=good,local_directional_stability=good))
    def test_first_running_blind_chain_holds_other_origins(self):
        self.add('screen_A');self.add('global_A','RUNNING');driver.advance(self.s)
        self.assertEqual(self.submitted,[])
    def test_complete_first_chain_releases_all_three_without_throttle(self):
        self.chain('A');driver.advance(self.s)
        self.assertEqual(self.submitted,['screen_B','screen_C','screen_D'])
    def test_failed_parent_blocks_new_dispatch(self):
        self.add('screen_A','FAILED')
        with self.assertRaisesRegex(ValueError,'HOLD'):driver.advance(self.s)
        self.assertEqual(self.submitted,[])
    def test_positive_control_does_not_count_as_blind_recovery(self):
        for origin in 'ABCD':self.chain(origin,good=origin=='A')
        driver.advance(self.s);self.assertEqual(self.submitted,['screen_E','screen_F'])
    def test_at_most_two_reserve_origins_then_unresolved(self):
        for origin in 'ABCDEF':self.chain(origin,good=False)
        driver.advance(self.s)
        self.assertEqual(self.submitted,[]);self.assertEqual(self.finished,['optimizer_unresolved'])
    def test_two_independent_passes_close_without_unneeded_reserve(self):
        for origin in 'ABCD':self.chain(origin,good=origin in 'AB')
        driver.advance(self.s)
        self.assertEqual(self.submitted,[]);self.assertEqual(self.finished,['conditional_recovery_pass'])
    def test_reviewed_failed_positive_is_replaced_by_actual_retry_task(self):
        self.rows=[r for r in self.rows if r['mode']!='positive_check'];self.add('positive_check','FAILED');self.add('positive_check',task='positive_check_retry1')
        w.table(self.s/'tables/reviewed_failures.tsv',[dict(task='positive_check')],['task'])
        w.typed(self.s/'checks/positive_check_retry1/result.tsv',dict(prediction_recovery=True,local_directional_stability=True))
        self.add('screen_A');self.add('global_A','RUNNING');driver.advance(self.s);self.assertEqual(self.submitted,[])
if __name__=='__main__':unittest.main()
