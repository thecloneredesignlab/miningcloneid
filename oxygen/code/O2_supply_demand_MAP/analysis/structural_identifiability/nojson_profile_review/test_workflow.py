#!/usr/bin/env python3
"""Metadata/controller tests; no model, Slurm submission or forbidden format."""
import multiprocessing
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch
import workflow as W
L=W.L

def job(tid='x',role='local',phase='primary',state='RUNNING',actual='',reserved='6'):
    return dict(slot=tid,task=tid,mode=role,role=role,phase=phase,job_id='123',token='owned_'+tid,state=state,
                exit_code='0:0' if state=='COMPLETED' else '',cpus='1',wall_seconds='21600',reserved_ch=reserved,actual_ch=actual,submitted_utc=L.now())

def increment(path):
    s=L.State(path)
    for _ in range(8):
        with s.locked():
            d=s.load();s.event(d,'concurrent','test');s.save(d)

class Tests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.p=Path(self.tmp.name)
        for sub in ['state','tasks','logs','submissions','tables','fits','manifests','reserve']: (self.p/sub).mkdir()
        self.s=L.State(self.p)
        self.cases=[dict(point=p,target='p_misseg' if p.startswith('pm') else 'buffer_beta',value=v) for p,v in [('pm_00',0),('bb_00',0),('pm_05',.01),('bb_07',1.6),('pm_06',.011),('bb_08',2.68),('pm_09',.4)]]
        L.write_table(self.p/'cases.tsv',self.cases,['point','target','value'])
    def tearDown(self):self.tmp.cleanup()
    def test_atomic_generation_and_corruption(self):
        d=self.s.load();d['reservations'].append(dict(task='paired',role='local',phase='primary',reserved_ch='6'));self.s.save(d)
        (self.p/'state/.partial').mkdir();self.assertEqual(self.s.load()['reservations'],d['reservations'])
        ptr=(self.p/'state/CURRENT').read_text().strip()
        with (self.p/'state'/ptr/'reservations.tsv').open('a') as f:f.write('corrupt')
        with self.assertRaises(ValueError):self.s.load()
    def test_concurrent_reducers(self):
        ps=[multiprocessing.Process(target=increment,args=(self.p,)) for _ in range(3)]
        for p in ps:p.start()
        for p in ps:p.join();self.assertEqual(p.exitcode,0)
        self.assertEqual(len(self.s.load()['events']),24)
    def test_audit_serializes_with_generation_publication(self):
        (self.p/'implementation').mkdir()
        L.write_table(self.p/'manifests/code_manifest.tsv',[],['path','sha256'])
        L.write_table(self.p/'manifests/control_manifest.tsv',[],['path','sha256'])
        writer=multiprocessing.Process(target=increment,args=(self.p,));writer.start()
        for _ in range(20):W.audit_locked(self.p)
        writer.join();self.assertEqual(writer.exitcode,0);self.assertEqual(len(self.s.load()['events']),8)
    def test_total_nonfit_and_paired_limits(self):
        d=self.s.load();d['jobs']=[job(reserved='3999')]
        with self.assertRaises(ValueError):L.budget_check(d,'local',1,21600)
        d['jobs']=[job(role='check',reserved='120')]
        with self.assertRaises(ValueError):L.budget_check(d,'check',1,1800)
        d['jobs']=[job(reserved='3900')]
        with self.assertRaises(ValueError):W.pair_reservation(self.s,d,dict(task_id='a_global',phase='primary'))
        d['jobs']=[];W.pair_reservation(self.s,d,dict(task_id='a_global',phase='primary'))
        self.assertEqual(L.committed(d),6);self.assertEqual(W.counts(d),(0,1))
    def test_pair_replaced_atomically_by_local_intent(self):
        d=self.s.load();d['jobs']=[job(str(i),state='COMPLETED',actual='0') for i in range(69)]
        d['reservations']=[dict(task='paired_local',role='local',phase='primary',reserved_ch='6')]
        with patch.object(L.subprocess,'run',return_value=subprocess.CompletedProcess([],0,'456\n','')):
            W.submit(self.s,d,self.p,'paired_local','local','local','primary')
        self.assertEqual(W.counts(d)[1],70);self.assertEqual(self.s.load()['reservations'],[])
        with self.assertRaises(ValueError):W.submit(self.s,d,self.p,'extra','local','local','primary')
    def test_failed_screens_consume_unique_origins(self):
        d=self.s.load();d['jobs']=[job(str(i),role='screen',state='FAILED',actual='0',reserved='1') for i in range(36)]
        self.assertEqual(W.counts(d)[0],36)
        with self.assertRaises(ValueError):W.submit(self.s,d,self.p,'extra','screen','screen','primary')
        d['jobs'].append(job('0',role='global',state='COMPLETED',actual='0',reserved='96'))
        self.assertEqual(W.counts(d)[0],36)
    def test_round_limits_and_no_rollover(self):
        d=self.s.load();d['jobs']=[job(str(i),role='screen',phase='reserve_1',state='FAILED',actual='0') for i in range(4)]
        with self.assertRaises(ValueError):W.submit(self.s,d,self.p,'x','screen','screen','reserve_1')
        d['jobs']=[job(str(i),phase='reserve_2',state='COMPLETED',actual='0') for i in range(8)]
        with self.assertRaises(ValueError):W.submit(self.s,d,self.p,'x','local','local','reserve_2')
        for n in (1,2):(self.p/'reserve'/str(n)).mkdir();(self.p/'reserve'/str(n)/'CLOSED').write_text('closed')
        self.assertIsNone(W.round_number(self.p))
    def test_ambiguous_submission_is_not_retried(self):
        d=self.s.load()
        with patch.object(L.subprocess,'run',side_effect=subprocess.TimeoutExpired('sbatch',90)) as call:
            with self.assertRaises(subprocess.TimeoutExpired):W.submit(self.s,d,self.p,'x','local','local','primary')
            W.submit(self.s,d,self.p,'x','local','local','primary');self.assertEqual(call.call_count,1)
            self.assertIn('--no-requeue',call.call_args[0][0]);self.assertIn('--qos=xxlarge',call.call_args[0][0])
        responses=[subprocess.CompletedProcess([],0,'456|'+d['jobs'][0]['token']+'|\n',''),subprocess.CompletedProcess([],0,'456|COMPLETED|0:0|1|20|20|\n',''),subprocess.CompletedProcess([],0,'456|COMPLETED|0:0|1|20|20|100M|\n','')]
        with patch.object(L.subprocess,'run',side_effect=responses):L.refresh(self.s,d)
        self.assertAlmostEqual(L.committed(d),20/3600)
    def test_task_collision_and_typed_roundtrip(self):
        v=dict(zero=0.,number=1.2345678901234567,flag=False,missing=None,label='English\tquote"')
        L.write_kv(self.p/'typed.tsv',v);self.assertEqual(L.read_kv(self.p/'typed.tsv'),v)
        t={k:None for k in L.TASK_FIELDS};t.update(task_id='a',mode='local');d=self.s.load()
        W.add(self.s,d,t);W.add(self.s,d,t);t['mode']='global'
        with self.assertRaises(ValueError):W.add(self.s,d,t)
    def test_partial_stale_and_skipped_acceptance(self):
        d=self.s.load();p=self.p/'fits/a';p.mkdir();(p/'accepted_fit.rds').write_text('opaque')
        self.assertIsNone(L.accepted(d,self.p,'a'))
        (p/'COMPLETE').write_text('done');L.write_kv(p/'result.tsv',dict(objective=3.8,skipped=False))
        (p/'payload_manifest.tsv').write_text('manifest')
        def receipt(fit):L.write_kv(p/'artifact_validation.tsv',dict(pass_gate=True,**{'pass':True},result_sha256=L.sha(p/'result.tsv'),accepted_fit_sha256=fit,payload_manifest_sha256=L.sha(p/'payload_manifest.tsv')))
        receipt(L.sha(p/'accepted_fit.rds'));d['jobs']=[dict(task='a',mode='check',state='COMPLETED',exit_code='0:0')]
        self.assertEqual(L.accepted(d,self.p,'a')['objective'],3.8)
        (p/'accepted_fit.rds').write_text('mutated')
        with self.assertRaises(ValueError):L.accepted(d,self.p,'a')
        L.write_kv(p/'result.tsv',dict(skipped=True));receipt(None)
        with self.assertRaises(ValueError):L.accepted(d,self.p,'a')
        (p/'TRANSFER_SKIPPED').write_text('all probes infeasible');self.assertTrue(L.accepted(d,self.p,'a')['skipped'])
    def candidate_rows(self):
        def row(tid,point,q,origin='',family='',mode='import',phase='historical_point',start=None):
            return dict(task_id=tid,point=point,objective=q,origin=origin or tid,family=family,mode=mode,phase=phase,alternative=phase in {'guard','repeat'},local_stable=True,source_fit=tid+'.rds',start_objective=q if start is None else start)
        rows=[row('reference','base',3),row('base_A','base',3.0001,'base_A','A'),row('base_B','base',3.0002,'base_B','B'),row('guard','base',3,mode='local',phase='guard',start=3.000001)]
        for c in self.cases:
            p=c['point'];rows.extend([row(p+'_old',p,4),row(p+'_A_local',p,4.00001,p+'_A','A','local','primary'),row(p+'_B_local',p,4.00002,p+'_B','B','local','primary'),row(p+'_repeat',p,4,mode='local',phase='repeat',start=4.000001)])
        return rows
    def test_cross_family_complete_chain_and_final_repeat(self):
        rows=self.candidate_rows()
        with patch.object(W,'candidates',return_value=rows),patch.object(W,'accepted',return_value={'objective':4}):
            bg,pg,cs=W.gates(self.s.load(),self.p,self.p);self.assertTrue(bg['pass']);self.assertTrue(all(r['pass_gate'] for r in pg))
            for r in rows:
                if r['point']=='pm_00' and r['family']=='B':r['origin']='pm_00_A';r['family']='A'
            self.assertFalse(W.gates(self.s.load(),self.p,self.p)[1][0]['independent_pass'])
    def test_better_warm_candidate_resets_point_recovery(self):
        rows=self.candidate_rows();rows[4]['objective']=3.8
        with patch.object(W,'candidates',return_value=rows),patch.object(W,'accepted',return_value={'objective':4}):
            pg=W.gates(self.s.load(),self.p,self.p)[1];self.assertFalse(pg[0]['pass_gate']);self.assertEqual(pg[0]['objective'],3.8)
    def test_global_without_accepted_parent_cannot_count(self):
        with patch.object(W,'candidates',return_value=self.candidate_rows()),patch.object(W,'accepted',return_value=None):
            self.assertTrue(all(not r['independent_pass'] for r in W.gates(self.s.load(),self.p,self.p)[1]))
    def test_ordinary_polish_is_not_final_repeat(self):
        rows=self.candidate_rows()
        for r in rows:
            if r['point']=='pm_00' and r['phase']=='repeat':r.update(phase='reserve_1',alternative=False)
        with patch.object(W,'candidates',return_value=rows),patch.object(W,'accepted',return_value={'objective':4}):
            self.assertFalse(W.gates(self.s.load(),self.p,self.p)[1][0]['repeat_pass'])
    def test_final_nonfit_closure_allowance_is_protected(self):
        d=self.s.load();d['jobs']=[job(role='check',reserved='119.5')]
        with self.assertRaises(ValueError):L.budget_check(d,'check',1,1800)
        self.assertEqual(L.budget_check(d,'report',1,1800),.5)
    def test_no_early_diagnostic_snapshot_while_fit_is_live(self):
        t={k:None for k in L.TASK_FIELDS};t.update(task_id='live',mode='global',phase='primary')
        d=self.s.load();W.add(self.s,d,t);d['jobs']=[job('live',role='global')];self.s.save(d)
        with patch.object(W,'drive',return_value='complete_unresolved_points'),patch.object(W,'emit'),patch.object(W,'ready',return_value=False):W.controller(self.p,self.p,True)
        self.assertFalse((self.p/'manifests/diagnostic_sources.tsv').exists())
    def test_format_guard_and_actual_resource_roles(self):
        L.write_table(self.p/'bad.tsv',[dict(value='[0, invalid transport')],['value'])
        with self.assertRaises(ValueError):L.audit(self.p)
        self.assertEqual(L.resources('global'),(8,43200,16));self.assertEqual(L.resources('local'),(1,21600,8))

if __name__=='__main__':unittest.main(verbosity=2)
