import unittest
import collect as c

def row(path,size=1):return dict(path=path,bytes=str(size),sha256='fixture')
class CollectionTests(unittest.TestCase):
    def test_full_accepted_states_and_actual_starts_are_kept(self):
        for path in ('fits/global_A/accepted_candidate.rds','checks/check_A/verified_candidate.rds','fits/global_A/actual_initial_population.rds','truth/truth_parameters.rds','reports/science/parameter_recovery.tsv'):
            self.assertTrue(c.eligible(row(path)))
    def test_large_redundant_histories_stay_canonical_on_HPC(self):
        for path in ('fits/global_A/improvements/123_0001234.rds','fits/global_A/best_123.rds','runtime/global_A/cache.so'):
            self.assertFalse(c.eligible(row(path)))
    def test_unsafe_or_forbidden_paths_fail_closed(self):
        for path in ('../outside','/absolute','state.json','view.js'):
            with self.assertRaises(ValueError):c.eligible(row(path))
    def test_disk_reserve_is_admitted_before_transfer(self):
        with self.assertRaisesRegex(ValueError,'reserve'):c.admitted([row('truth/data.rds',20)],c.RESERVE+19)
        selected,total=c.admitted([row('truth/data.rds',20)],c.RESERVE+20)
        self.assertEqual((len(selected),total),(1,20))
if __name__=='__main__':unittest.main()
