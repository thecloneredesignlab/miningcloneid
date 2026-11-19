"""Summarize sealed D19 outputs and independent audits without model calls."""
import csv
from pathlib import Path
import sys

from collect_closeout import ORIGINS, PREDICTION_TERMINALS, table, truth

def write_tsv(path, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)

def closeout(root, out):
    budget, = table(root/'accounting/final_budget.tsv')
    terminal = (root/'controller/TERMINAL').read_text().strip()
    assert budget['status']==terminal
    lines = ['# D19 deterministic blind recovery closeout', '',
             f'**{terminal}**', '',
             f"Settled allocation use was {float(budget['settled_core_hours']):.6f} of 120 core-hours. "
             f"The {budget['owned_jobs']} owned allocations had zero live jobs at reconciliation. "
             f"Charged forwards totaled {int(budget['attempts']):,} of 336,400. "
             f"Reserve rounds used: {budget['scientific_reserve_rounds']} of one. "
             f"Technical retries used: {budget['technical_retries']} of one.", '']
    decision = dict(d19_status=terminal, prediction_gate='WITHHELD', lowest_D_local_quality='WITHHELD',
                    scientific_n20_review='HELD_D19_AUDIT', released_N20=0, released_profiles=0, released_tumor=0)
    retry_path = root/'accounting/technical_retry.tsv'
    if retry_path.is_file():
        retry, = table(retry_path)
        lines += [f"B's original second replay ({retry['previous_job']}) stopped when the clean container "
                  'environment omitted the Slurm allocation ID. The saved coordinates, selections and predictions '
                  'agreed exactly with the first replay. The single charged retry forwarded the native ID through '
                  'the container and reused the previous verified checkpoint. The original failed call and allocation '
                  'remain in the totals. The original release and scientific code remain sealed.', '']
    if terminal not in PREDICTION_TERMINALS:
        record, = table(root/'accounting/stop_record.tsv')
        lines += [record['detail'], '', 'Report scope is limited to the sealed records available at the stop.', '']
    else:
        quality = table(root/'reports/verified_candidates.tsv')
        verdicts = table(root/'reports/origin_verdicts.tsv')
        replays = table(out/'d19_replay_gate_audit.tsv')
        polls = table(out/'d19_local_poll_audit.tsv')
        lock, = table(root/'reports/selection_lock.tsv')
        assessment, = table(root/'verification/P/assessment_gate.tsv')
        if assessment['status']!='PASS' or int(assessment['full_forward_calls']) != 0:
            raise RuntimeError('Post-selection assessor receipt invalid')
        gates = {}
        for row in quality:
            origin, key = row['origin'], row['key']
            pair = [r for r in replays if r['origin']==origin and r['key']==key]
            if len(pair)!=2 or not all(truth(r['passes_all_four'])==truth(row['prediction_pass']) for r in pair):
                raise RuntimeError('Independent replay gates disagree: '+origin+'/'+key)
            if any(truth(r['passes_all_four']) for r in pair):
                gates[(origin,key)] = row

        outcomes = []
        for row in verdicts:
            origin = row['origin']
            candidates = [r for r in quality if r['origin']==origin]
            if not candidates:
                raise RuntimeError('Missing verified coordinates for origin '+origin)
            best = min(candidates,key=lambda r:float(r['D']))
            eligible = [r for (o,k),r in gates.items() if o==origin]
            strict = min(eligible,key=lambda r:float(r['D'])) if eligible else None
            assert bool(strict)==truth(row['prediction_pass'])
            if strict is not None:
                assert row['strict_key']==strict['key']
            points = [('lowest_verified_D',best)] + ([('lowest_D_strict_pass',strict)] if strict else [])
            for label, point in points:
                replay = next(r for r in replays if r['origin']==origin and r['key']==point['key'])
                outcomes.append(dict(origin=origin, vector_role=label, key=point['key'], best_found_key=row['best_found_key'],
                    D=point['D'], D_G=point['D_G'], D_K=point['D_K'], growth_RMS_over_SD=point['rms_growth'],
                    max_growth_over_SD=point['max_growth'], max_karyotype_TV=point['max_kary_tv'],
                    worst_karyotype_sample=replay['worst_sample'], independent_replays=2,
                    D_gate=replay['D_gate'], growth_RMS_gate=replay['growth_RMS_gate'],
                    growth_max_gate=replay['growth_max_gate'], karyotype_gate=replay['karyotype_gate'],
                    all_prediction_gates=point['prediction_pass'], local_poll_pass=point['local_pass']))
        write_tsv(out/'d19_origin_and_gate_audit.tsv',outcomes)
        eligible_origins = set(origin for origin,_ in gates)
        num = len(eligible_origins)
        assert int(lock['independent_prediction_origins'])==num
        f_used = any(r['origin']=='D_F' for r in verdicts)
        assert f_used==bool(int(budget['scientific_reserve_rounds']))
        assert f_used==(len(eligible_origins & set(ORIGINS[:3]))<2)
        lowest = min(verdicts,key=lambda r:float(r['best_found_D']))
        assert (lowest['origin'],lowest['best_found_key']) == (lock['best_found_origin'],lock['best_found_key'])
        lowest_poll = [r for r in polls if truth(r['is_global_lowest_found']) and truth(r['associated_final_poll'])]
        local = len(lowest_poll)==1 and truth(lowest_poll[0]['local_pass'])
        assert local==truth(lock['best_found_local_pass'])
        assert terminal == ('DET_OPTIMIZER_UNRESOLVED' if num < 2 else
           'D19_PREDICTION_RECOVERY_PASS' if local else 'D19_PREDICTION_RECOVERED_LOCAL_UNRESOLVED')
        decision.update(prediction_gate='PASS' if num>=2 else 'UNRESOLVED',
                        lowest_D_local_quality='PASS' if local else 'UNRESOLVED',
                        scientific_n20_review='ELIGIBLE_FOR_SEPARATE_REVIEW' if num>=2 and local else
                        'REVIEW_D19_LOCAL_QUALITY' if num>=2 else 'HELD_D19_UNRESOLVED')
        origin_list = ', '.join(sorted(eligible_origins)) or 'none'
        lines += [f'{num} independent origins passed all four criteria in two fresh replays each. Passing origins: {origin_list}.',
                  'The isolated generating-vector control passed D <= 1e-12 and receives zero blind-origin credit.', '',
                  '| Origin | Vector | D | Growth RMS / SD | Max growth / SD | Max sample TV | Strict prediction pass |',
                  '|---|---|---:|---:|---:|---:|---|']
        for r in outcomes:
            label = r['vector_role'].replace('_',' ')
            lines.append(f"| {r['origin']} | {label} | {float(r['D']):.8g} | {float(r['growth_RMS_over_SD']):.6g} | "
                         f"{float(r['max_growth_over_SD']):.6g} | {float(r['max_karyotype_TV']):.6g} | {r['all_prediction_gates']} |")
        lines += ['',f"The exact lowest D was {float(lowest['best_found_D']):.12g} from {lowest['origin']}. "
                     f"Vector key: `{lowest['best_found_key']}`.",
                  f"Independent replay of this exact vector: {lock['best_found_replay_verified']}. "
                  f"Local poll qualification: {local}."]
        if lowest_poll:
            poll = lowest_poll[0]
            lines.append(f"The anchor had {poll['evaluated_directions']}/{poll['expected_directions']} registered directions "
                         f"with maximum improvement {float(poll['max_D_improvement']):.8g}.")
        else:
            lines.append('The actual lowest D has no complete registered exact-anchor poll certificate.')
        lines += ['', 'Per-sample TV and exact replay/poll audits accompany the origin table.',
                      'The post-selection assessment preserves all 20 roles, including the fixed generating growth scale, '
                      'and reports raw widths, kernels, root distributions and mechanism-response error separately.', '']
    lines += ['This is a conditional deterministic culture test on one fixed seed144 case. '
              'The result scope is prediction recovery and finite numerical optimization under the declared passage rules.',
              'The actual sparse-design structural question and repeated noisy recovery remain open.',
              'N20, profile and tumor allocations released by this audit: zero.', '']
    (out/'d19_closeout_review.md').write_text('\n'.join(lines))
    write_tsv(out/'N20_release_decision.tsv',[decision])
    write_tsv(out/'d19_recovery_scope.tsv',[
        dict(topic='deterministic_growth_and_karyotype',status=decision['prediction_gate'],
             scope='Frozen seed144 conditioned observation design'),
        dict(topic='finite_exact_vector_local_optimization',status=decision['lowest_D_local_quality'],
             scope='Lowest found D, registered directions and meshes only'),
        dict(topic='actual_sparse_passage_structural_identifiability',status='OPEN',
             scope='Earlier ideal structural witness results retain their original scope'),
        dict(topic='noisy_repeated_recovery',status='OPEN',scope='Separate noisy validation and allocation')])

if __name__=='__main__':
    closeout(Path(sys.argv[1]),Path(sys.argv[2]))
