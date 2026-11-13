"""Mocked controller fixtures only; never invokes Slurm or the model."""
import sys,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
import workflow as w
import profile_windows as p

class ProfileDispatch(unittest.TestCase):
 def setUp(self):
  self.tmp=tempfile.TemporaryDirectory(prefix='ivc_profile_metadata_fixture_');self.s=Path(self.tmp.name)
  for name in ('tables','manifests','profile_epochs/e0001'):(self.s/name).mkdir(parents=True)
  (self.s/'manifests/input_manifest.tsv').write_text('path\tsha256\n')
  (self.s/'profile_epochs/e0001/manifest.tsv').write_text('path\tsha256\n')
  for policy in p.POLICIES:(self.s/'profile_epochs/e0001'/('baseline_'+policy+'.rds')).write_bytes(b'LABELED_METADATA_FIXTURE_NOT_RDS')
  w.typed(self.s/'tables/profile_epoch.tsv',dict(epoch='e0001',manifest='profile_epochs/e0001/manifest.tsv',manifest_sha256=w.sha(self.s/'profile_epochs/e0001/manifest.tsv'),input_manifest_sha256=w.sha(self.s/'manifests/input_manifest.tsv')))
  self.set_gate(False)
  w.table(self.s/'tables/profile_anchor_gates.tsv',[dict(policy=policy,index=index,pass_value='false') for policy in p.POLICIES for index in range(1,8)],['policy','index','pass_value'])
  text=(self.s/'tables/profile_anchor_gates.tsv').read_text().replace('pass_value','pass');(self.s/'tables/profile_anchor_gates.tsv').write_text(text)
  self.rows=[];self.sent=[]
 def tearDown(self):self.tmp.cleanup()
 def set_gate(self,invalid):w.typed(self.s/'tables/profile_gate.tsv',{'pass':False,'free_baseline_invalidated':invalid,'evidence_manifest':'tables/fixture.tsv'})
 def submit(self,study,task,**options):
  self.assertFalse(any(r['task']==task for r in self.rows))
  j=dict(options,task=task,state='SUBMITTED',exit_code='',source_task=options.get('source_task',''))
  self.rows.append(j);self.sent.append(j)
 def run_pass(self,baseline_ok=True):
  with patch.object(sys,'argv',['profile_windows.py',str(self.s)]),patch.object(w,'snapshot',side_effect=lambda s:list(self.rows)),patch.object(w,'gate',side_effect=None if baseline_ok else ValueError('baseline not accepted')),patch.object(w,'verify_manifest'),patch.object(w,'submit',side_effect=self.submit):p.main()
 def test_ready_design_and_idempotency(self):
  self.run_pass();self.assertEqual(len(self.sent),18)
  self.assertEqual(sum(j['mode']=='global' for j in self.sent),14)
  self.assertEqual({j['role'] for j in self.sent if j['mode']=='local'},{'density_ascending','density_descending'})
  self.assertTrue(all(j['kind']=='profile' and j['root_family']=='profile' for j in self.sent))
  self.run_pass();self.assertEqual(len(self.sent),18)
 def test_baseline_gate_blocks_all(self):
  with self.assertRaisesRegex(ValueError,'baseline not accepted'):self.run_pass(False)
  self.assertFalse(self.sent)
 def test_counterexample_blocks_all(self):
  self.set_gate(True)
  with self.assertRaisesRegex(ValueError,'invalidated'):self.run_pass()
  self.assertFalse(self.sent)
 def test_no_iteration_or_attempt_cap(self):
  self.assertEqual(p.next_id('fixture_w9999'),'fixture_w10000')
  self.assertEqual(p.next_local_decision({'termination':'time_limit'},{'directional_improvement':'0'},False),'checkpoint')
  self.assertEqual(p.next_local_decision({'termination':'call_limit'},{'directional_improvement':'.003'},False),'endpoint')
  self.assertEqual(p.next_local_decision({'termination':'poll_stationary'},{'directional_improvement':'0'},False),'stationary')
 def test_cross_branch_and_nonadjacent_winners_rejected(self):
  parent=dict(task='a',kind='profile',mode='local',role='density_ascending',policy='production',root_id='profile_e0001_production_ascending',density_index='1')
  options=dict(mode='local',role='density_ascending',policy='production',root_id=parent['root_id'],density_index=3,selector='task',source_task='a')
  with self.assertRaisesRegex(ValueError,'adjacent'):p.verify_branch_source(self.s,'e0001','next_w0001',options,{'a':parent})
  options['density_index']=2;p.verify_branch_source(self.s,'e0001','next_w0001',options,{'a':parent})
  parent['root_id']='profile_e0001_production_descending'
  with self.assertRaisesRegex(ValueError,'Cross-branch'):p.verify_branch_source(self.s,'e0001','next_w0001',options,{'a':parent})
 def test_material_probe_returns_to_polish_before_next_anchor(self):
  self.run_pass();j=next(r for r in self.rows if r['policy']=='production' and r['role']=='density_ascending');j.update(state='COMPLETED',exit_code='0:0')
  task=j['task'];check=task+'_check';self.rows.append(dict(j,task=check,mode='profile_fresh',role='density_check',source_task=task))
  (self.s/'fits'/task).mkdir(parents=True);(self.s/'tests'/check).mkdir(parents=True)
  w.typed(self.s/'fits'/task/'result.tsv',{'termination':'poll_stationary'})
  w.typed(self.s/'tests'/check/'result.tsv',{'directional_improvement':.003})
  w.table(self.s/'tables/profile_branch_history.tsv',[dict(policy='production',index=1,branch='ascending',objective=10,source_task=check,termination='poll_stationary',directional_gain=.003)],['policy','index','branch','objective','source_task','termination','directional_gain'])
  self.run_pass();new=self.sent[18:]
  self.assertEqual(len(new),1);self.assertEqual(new[0]['task'],p.next_id(task));self.assertEqual(new[0]['source_task'],check)
  self.assertFalse(any(r['task']=='p_e0001_production_r02_ascending_w0001' for r in new))
if __name__=='__main__':unittest.main()
