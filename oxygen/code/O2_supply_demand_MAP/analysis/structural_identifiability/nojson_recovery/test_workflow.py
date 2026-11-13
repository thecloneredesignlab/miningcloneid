#!/usr/bin/env python3
import importlib.util
import multiprocessing
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import subprocess
import sys
spec=importlib.util.spec_from_file_location('workflow',Path(__file__).with_name('workflow.py'))
w=importlib.util.module_from_spec(spec); spec.loader.exec_module(w)

def increment(path):
    s=w.State(path)
    for _ in range(8):
        with s.locked():
            d=s.load(); s.event(d,'concurrent','record'); s.save(d)

class Tests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(); self.p=Path(self.tmp.name)
        for sub in ['state','tasks','logs','submissions','tables','fits']: (self.p/sub).mkdir()
        self.s=w.State(self.p)
    def tearDown(self): self.tmp.cleanup()
    def test_atomic_partial_and_corruption(self):
        d=self.s.load(); self.s.event(d,'first','ready'); self.s.save(d)
        (self.p/'state/.partial').mkdir()
        self.assertEqual(len(self.s.load()['events']),1)
        ptr=(self.p/'state/CURRENT').read_text().strip()
        with (self.p/'state'/ptr/'jobs.tsv').open('a') as f: f.write('partial')
        with self.assertRaises(ValueError): self.s.load()
    def test_concurrent_reducers(self):
        ps=[multiprocessing.Process(target=increment,args=(self.p,)) for _ in range(3)]
        for p in ps: p.start()
        for p in ps: p.join(); self.assertEqual(p.exitcode,0)
        self.assertEqual(len(self.s.load()['events']),24)
    def test_budget_reservations_and_caps(self):
        d=self.s.load()
        for i in range(15):
            d['jobs'].append(dict(role='global',state='RUNNING',reserved_ch='96',actual_ch=''))
        with self.assertRaises(ValueError): w.budget_check(d,'local',1,21600)
        d['jobs']=d['jobs'][:1]; d['jobs'][0].update(state='COMPLETED',actual_ch='1')
        self.assertEqual(w.committed(d),1)
        d['jobs']=[dict(role='check',state='RUNNING',reserved_ch='120',actual_ch='')]
        with self.assertRaises(ValueError): w.budget_check(d,'check',1,1800)
    def test_submission_intent_ambiguity_and_idempotence(self):
        d=self.s.load()
        with patch.object(w.subprocess,'run',side_effect=subprocess.TimeoutExpired('sbatch',90)) as call:
            with self.assertRaises(subprocess.TimeoutExpired): w.submit(self.s,d,self.p,'x','local','local','diagnostic')
            self.assertEqual(self.s.load()['jobs'][0]['state'],'AMBIGUOUS')
            w.submit(self.s,d,self.p,'x','local','local','diagnostic')
            self.assertEqual(call.call_count,1)
            self.assertIn("--no-requeue",call.call_args[0][0])
        responses=[subprocess.CompletedProcess([],0,'123|'+d['jobs'][0]['token']+'|\n',''),subprocess.CompletedProcess([],0,'123|COMPLETED|0:0|1|20|20|\n',''),subprocess.CompletedProcess([],0,'123|COMPLETED|0:0|1|20|20||\n123.batch|COMPLETED|0:0|1|20|20|100M|\n','')]
        with patch.object(w.subprocess,'run',side_effect=responses): w.refresh(self.s,d)
        self.assertEqual(d['jobs'][0]['job_id'],'123'); self.assertAlmostEqual(float(d['jobs'][0]['actual_ch']),20/3600)
    def test_unexpected_fit_failure_holds_dependents(self):
        d=self.s.load()
        d['jobs']=[dict(slot='x',task='failed',mode='local',role='local',phase='primary',job_id='456',token='owned',state='RUNNING',exit_code='',cpus='1',wall_seconds='21600',reserved_ch='6',actual_ch='',submitted_utc=w.now())]
        result=subprocess.CompletedProcess([],0,'456|FAILED|1:0|1|10|10|\n','')
        with patch.object(w.subprocess,'run',return_value=result):
            with self.assertRaisesRegex(ValueError,'dependent work held'): w.refresh(self.s,d)
        self.assertEqual(self.s.load()['jobs'][0]['state'],'FAILED')
    def test_real_allocation_policy(self):
        self.assertEqual(w.resources('global'),(8,43200,16))
        self.assertEqual(w.resources('local'),(1,21600,8))
    def test_accepted_payload_mutation(self):
        d=self.s.load(); target=self.p/'fits/a'; target.mkdir()
        w.atomic_text(target/'accepted_fit.rds','opaque test payload')
        w.write_kv(target/'result.tsv',{'objective':3.8})
        w.atomic_text(target/'payload_manifest.tsv','test payload manifest')
        w.atomic_text(target/'COMPLETE','complete')
        w.write_kv(target/'artifact_validation.tsv',{'pass':True,'result_sha256':w.sha(target/'result.tsv'),'accepted_fit_sha256':w.sha(target/'accepted_fit.rds'),'payload_manifest_sha256':w.sha(target/'payload_manifest.tsv')})
        d['jobs']=[dict(task='a',mode='check',state='COMPLETED',exit_code='0:0')]
        self.assertEqual(w.accepted(d,self.p,'a')['objective'],3.8)
        w.atomic_text(target/'accepted_fit.rds','changed payload')
        with self.assertRaises(ValueError): w.accepted(d,self.p,'a')
    def test_partial_not_accepted(self):
        d=self.s.load(); (self.p/'fits/a').mkdir(); (self.p/'fits/a/accepted_fit.rds').touch()
        self.assertIsNone(w.accepted(d,self.p,'a'))
    def test_stale_snapshot(self):
        (self.p/'implementation').mkdir(); (self.p/'manifests').mkdir(); f=self.p/'implementation/a.R'; f.write_text('x <- 1\n')
        w.write_table(self.p/'manifests/code_manifest.tsv',[dict(path='a.R',sha256=w.sha(f))],['path','sha256'])
        w.write_table(self.p/'manifests/control_manifest.tsv',[],['path','sha256'])
        w.verify_code(self.p); f.write_text('x <- 2\n')
        with self.assertRaises(ValueError): w.verify_code(self.p)
    def test_format_guard(self):
        w.write_table(self.p/'good.tsv',[dict(value='ordinary text')],['value']); w.audit(self.p)
        w.write_table(self.p/'bad.tsv',[dict(value='[0, invalid serialized transport')],['value'])
        with self.assertRaises(ValueError): w.audit(self.p)
    def test_task_identity_and_typed_roundtrip(self):
        v=dict(zero=0.,integer=3,flag=False,missing=None,label='English\tquote"')
        w.write_kv(self.p/'typed.tsv',v); self.assertEqual(w.read_kv(self.p/'typed.tsv'),v)
        d=self.s.load(); t=w.task_record('a','local','diagnostic','old',source='old.rds')
        w.add_task(self.s,d,t); w.add_task(self.s,d,t); self.assertEqual(len(self.s.load()['tasks']),1)

if __name__=='__main__': unittest.main(verbosity=2)
