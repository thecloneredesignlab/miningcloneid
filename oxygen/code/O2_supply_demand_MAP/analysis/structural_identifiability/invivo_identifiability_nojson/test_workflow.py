#!/usr/bin/env python3
"""Control-plane tests; forbidden-format examples use mocked paths only."""
import tempfile
import unittest
from pathlib import Path
from unittest import mock
import workflow as w

class WorkflowTests(unittest.TestCase):
    def test_coherent_generation_and_corruption(self):
        with tempfile.TemporaryDirectory() as d:
            s=w.State(Path(d))
            with s.locked(): s.save([]); self.assertEqual(s.load(),[])
            current=(Path(d)/'state/CURRENT').read_text().strip()
            (Path(d)/'state'/current/'jobs.tsv').write_text('corrupted')
            with self.assertRaises(ValueError): s.load()
    def test_ambiguous_submission_and_reconciliation(self):
        rows=[dict(job_id='',comment='owned',state='INTENT')]
        with self.assertRaises(ValueError): w.resolve_matches(rows,[('12','other')])
        with self.assertRaises(ValueError): w.resolve_matches(rows,[('12','owned'),('13','owned')])
        w.resolve_matches(rows,[('12','owned'),('12','owned'),('14','other')])
        self.assertEqual(rows[0]['job_id'],'12')
        w.resolve_matches(rows,[('13','owned')]); self.assertEqual(rows[0]['job_id'],'12')
    def test_format_guard_without_writing_forbidden_fixture(self):
        item=mock.Mock(); item.is_file.return_value=True; item.suffix='.json'
        with mock.patch.object(Path,'rglob',return_value=[item]):
            with self.assertRaises(ValueError): w.audit_formats(Path('/unused'))
    def test_attempt_caps_and_reservations(self):
        b=[dict(role='global',attempts='2',cores='8',wall_seconds='28800',cap_core_hours='128')]
        self.assertEqual(w.check_budget([],b,'global')[1],64)
        job=dict(role='global',state='SUBMITTED',reserved_ch='64',actual_ch='')
        with self.assertRaises(ValueError): w.check_budget([job,job],b,'global')
        other=dict(role='check',state='SUBMITTED',reserved_ch='150',actual_ch='')
        with self.assertRaises(ValueError): w.check_budget([other],b,'global')
        failed=dict(role='global',state='FAILED',reserved_ch='64',actual_ch='1')
        self.assertEqual(w.charge(failed),1)
        with self.assertRaises(ValueError): w.check_budget([failed,failed],b,'global')
    def test_incomplete_and_stale_payload_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)
            with self.assertRaises(ValueError): w.verify_payload(p)
            (p/'record.tsv').write_text('value\n1\n')
            w.table(p/'payload_manifest.tsv',[dict(path='record.tsv',sha256=w.sha(p/'record.tsv'))],['path','sha256'])
            (p/'COMPLETE').write_text(w.sha(p/'payload_manifest.tsv')+'\n')
            w.verify_payload(p)
            (p/'record.tsv').write_text('value\n2\n')
            with self.assertRaises(ValueError): w.verify_payload(p)

if __name__=='__main__': unittest.main()
