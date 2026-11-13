#!/usr/bin/env python3
"""Exercise the full scientific dependency graph with an analytic fake fitter."""
import importlib.util
import tempfile
import unittest
import multiprocessing
import contextlib
import io
import sys
import time
from pathlib import Path

spec = importlib.util.spec_from_file_location("pilot", Path(__file__).with_name("culture_profile_pilot.py"))
pilot = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pilot)


class FakePilot(pilot.Pilot):
    def verify(self):
        pass

    def refresh(self):
        pass

    def ensure_controller(self, task):
        pass

    def check_saved_objects(self, task):
        path = self.study / "fits" / task["task_id"] / "artifact_validation_v05.json"
        p = pilot.read(self.study / "fits" / task["task_id"] / "validated_result.json")["params"]
        pilot.write(path, {"pass": True, "analytic_oracle_only": True, "numeric_params_exact": {n: format(v, '.17g') for n, v in p.items()}})
        return pilot.sha(path)

    def submit_command(self, args, directory):
        return str(9000 + self.state["next_ordinal"])

    def finish(self, task):
        starts = pilot.read(self.study / "tasks" / task["task_id"] / "starts.json")
        params = self.project_params(next(iter(starts.values())), task["fixed"])
        if task["stage"] == "baseline":
            q = 10. if task["reserve"] else 10. + .0001 * task["ordinal"]
        else:
            target, value = next(iter(task["fixed"].items()))
            q = 10. + (value - (.1 if target == "p_misseg" else 1.))**2
        result = dict(task_id=task["task_id"], params=params, objective=q, native_objective=q,
                      winner_phase="analytic_oracle", fit_wall_seconds=.1, total_evaluator_invocations=1,
                      allocated_core_hours=.001, search_status="analytic_test_only")
        pilot.write(self.study / "fits" / task["task_id"] / "validated_result.json", result)
        task["status"] = "accepted"


class CounterPilot:
    def __init__(self, project, study):
        self.study = Path(study)
        self.state = pilot.read(self.study / "state.json")
        time.sleep(.01)

    def advance(self, initial=False):
        self.state["counter"] += 1
        pilot.write(self.study / "state.json", self.state)
        pilot.write(self.study / "reports/progress.json", self.state)


def run_concurrent_reducers(folder):
    pilot.Pilot = CounterPilot
    for _ in range(5):
        sys.argv = ["pilot", "advance", folder, folder]
        with contextlib.redirect_stdout(io.StringIO()):
            pilot.main()


class DependencyGraph(unittest.TestCase):
    def make_pilot(self, root):
        pilot.write(root / "study.json", dict(master_seed=100, production_objective=10.1,
                                              zero_scales=dict(p_misseg=.015, buffer_beta=1.),
                                              grid_refinement_score_threshold=.05))
        labels = ["seed144", "seed371", "seed10", "engineering_null"] + ["independent_" + str(i) for i in range(1, 9)]
        for i, name in enumerate(labels):
            params = {n: 1. for n in pilot.SYMBOLS}
            params.update(lam_max=.1 + .01 * i, p_misseg=.1, buffer_beta=1.)
            pilot.write(root / "start_bank" / (name + ".json"), params)
        return FakePilot(root, root)

    def test_full_graph_submits_every_ready_start_and_finishes_within_work_budget(self):
        with tempfile.TemporaryDirectory() as folder:
            p = self.make_pilot(Path(folder))
            p.advance(initial=True)
            self.assertEqual(len(p.state["tasks"]), 5)
            for task in p.state["tasks"].values():
                p.finish(task)
            p.advance()
            self.assertEqual(len(p.state["tasks"]), 6)
            for task in p.state["tasks"].values():
                if not p.terminal(task):
                    p.finish(task)
            p.advance()
            ready = [x for x in p.state["tasks"].values() if x["stage"] == "profile" and not p.terminal(x)]
            self.assertGreaterEqual(len(ready), 80)
            for _ in range(40):
                for task in list(p.state["tasks"].values()):
                    if not p.terminal(task):
                        p.finish(task)
                p.advance()
                if p.state["status"] == "pilot_fitting_finished":
                    break
            self.assertEqual(p.state["status"], "pilot_fitting_finished")
            initial = [x for x in p.state["points"].values() if x["kind"] == "initial"]
            self.assertEqual(len(initial), 20)
            for point in initial:
                main = [p.state["tasks"][t] for t in point["task_ids"] if p.state["tasks"][t]["role"].startswith("point_start_")]
                self.assertEqual(len(main), 5)
            self.assertLessEqual(len(p.state["tasks"]), 125)
            self.assertLessEqual(p.state["fit_core_hours_reserved"], 8640)
            self.assertLessEqual(sum(x["reserve"] for x in p.state["tasks"].values()), 20)

    def test_better_constrained_solution_triggers_unrestricted_repair(self):
        with tempfile.TemporaryDirectory() as folder:
            p = self.make_pilot(Path(folder))
            p.advance(initial=True)
            for task in p.state["tasks"].values():
                p.finish(task)
            p.advance()
            for task in p.state["tasks"].values():
                if not p.terminal(task):
                    p.finish(task)
            p.advance()
            task = next(x for x in p.state["tasks"].values() if x["stage"] == "profile")
            p.finish(task)
            path = p.study / "fits" / task["task_id"] / "validated_result.json"
            value = pilot.read(path); value["objective"] = 9.9; pilot.write(path, value)
            p.advance()
            self.assertEqual(p.state["status"], "baseline_repair_running")
            self.assertTrue(any(x["role"] == "baseline_repair" for x in p.state["tasks"].values()))

    def test_opposite_direction_returns_from_null_to_nearest_positive(self):
        with tempfile.TemporaryDirectory() as folder:
            p = self.make_pilot(Path(folder)); p.advance(initial=True)
            for _ in range(30):
                for task in list(p.state["tasks"].values()):
                    if not p.terminal(task): p.finish(task)
                p.advance()
                if len(p.state["points"]) == 20 and all(
                    len([tid for tid in point["task_ids"] if p.state["tasks"][tid]["role"].startswith("point_start_") and p.terminal(p.state["tasks"][tid])]) == 5
                    for point in p.state["points"].values()): break
            null_checks = [t for t in p.state["tasks"].values() if t["role"] in ("null_reverse", "null_independent")]
            self.assertEqual(len(null_checks), 4)
            self.assertFalse(any(t["role"] == "point_reverse" for t in p.state["tasks"].values()))
            for task in null_checks: p.finish(task)
            p.reserve_checks()
            returns = [t for t in p.state["tasks"].values() if t["role"] == "point_reverse"]
            self.assertEqual(len(returns), 2)
            for task in returns:
                target, value = next(iter(task["fixed"].items()))
                positives = [x["value"] for x in p.state["points"].values() if x["target"] == target and x["value"] > 0]
                self.assertEqual(value, min(positives))
                starts = pilot.read(p.study / "tasks" / task["task_id"] / "starts.json")
                self.assertTrue(next(iter(starts)).startswith("reverse_from_null_"))
                self.assertEqual(next(iter(starts.values()))[target], 0.)
                self.assertTrue(task["reserve"])
            count = len(p.state["tasks"])
            p.reserve_checks()
            self.assertEqual(len(p.state["tasks"]), count)
            self.assertLessEqual(sum(t["reserve"] for t in p.state["tasks"].values()), 20)

    def test_reverse_return_does_not_exceed_remaining_shared_budget(self):
        with tempfile.TemporaryDirectory() as folder:
            p = self.make_pilot(Path(folder)); p.advance(initial=True)
            for _ in range(30):
                for task in list(p.state["tasks"].values()):
                    if not p.terminal(task): p.finish(task)
                p.advance()
                if len(p.state["points"]) == 20 and all(
                    len([tid for tid in point["task_ids"] if p.state["tasks"][tid]["role"].startswith("point_start_") and p.terminal(p.state["tasks"][tid])]) == 5
                    for point in p.state["points"].values()): break
            for task in list(p.state["tasks"].values()):
                if not p.terminal(task): p.finish(task)
            # Reserve bookkeeping only: leave one slot for the two ready returns.
            used = sum(t["reserve"] for t in p.state["tasks"].values())
            for i in range(19-used):
                p.state["tasks"]["spent_"+str(i)] = {"status":"failed", "reserve":True, "role":"spent"}
            p.reserve_checks()
            self.assertEqual(sum(t["role"] == "point_reverse" for t in p.state["tasks"].values()), 1)
            self.assertEqual(p.remaining_reserve(), 0)

    def test_warm_strengthening_cannot_replace_independent_agreement(self):
        with tempfile.TemporaryDirectory() as folder:
            p = self.make_pilot(Path(folder))
            p.advance(initial=True)
            for i, task in enumerate(p.state["tasks"].values()):
                p.finish(task)
                path = p.study / "fits" / task["task_id"] / "validated_result.json"
                value = pilot.read(path); value["objective"] = 10. + i; pilot.write(path, value)
            p.advance()
            strong = next(x for x in p.state["tasks"].values() if not p.terminal(x))
            p.finish(strong)
            self.assertFalse(strong["agreement_eligible"])
            p.advance()
            repeated = [x for x in p.state["tasks"].values() if x["role"] == "baseline_independent_repeat"]
            self.assertEqual(len(repeated), 5)
            self.assertTrue(all(x["agreement_eligible"] for x in repeated))
            self.assertFalse(p.state["points"])
            for task in repeated:
                starts = pilot.read(p.study / "tasks" / task["task_id"] / "starts.json")
                self.assertEqual(len(starts), 1)
                self.assertTrue(next(iter(starts)).startswith("independent_"))

    def test_concurrent_reducers_do_not_lose_state_updates(self):
        with tempfile.TemporaryDirectory() as folder:
            pilot.write(Path(folder) / "state.json", {"counter": 0})
            context = multiprocessing.get_context("fork")
            jobs = [context.Process(target=run_concurrent_reducers, args=(folder,)) for _ in range(4)]
            for job in jobs:
                job.start()
            for job in jobs:
                job.join(timeout=10)
                self.assertEqual(job.exitcode, 0)
            self.assertEqual(pilot.read(Path(folder) / "state.json")["counter"], 20)

    def test_refinement_does_not_promise_two_five_start_points_with_six_reserves(self):
        with tempfile.TemporaryDirectory() as folder:
            p = self.make_pilot(Path(folder))
            for target in pilot.TARGETS:
                for i, value in enumerate([0., .1, .4] if target == "p_misseg" else [0., 1., 10.]):
                    pid = target + "_" + str(i)
                    p.state["points"][pid] = dict(point_id=pid, target=target, value=value, kind="initial", task_ids=[])
            p.remaining_reserve = lambda: 6
            p.baseline = lambda: {"task_id": "analytic_baseline"}
            p.accepted = lambda point_id: [{"objective": 12. if point_id.endswith("_1") else 10.}]
            p.submit_points = lambda: None
            p.refine_grid()
            added = [point for point in p.state["points"].values() if point["kind"] == "refinement"]
            self.assertEqual(len(added), 1)


if __name__ == "__main__":
    unittest.main(verbosity=2)
