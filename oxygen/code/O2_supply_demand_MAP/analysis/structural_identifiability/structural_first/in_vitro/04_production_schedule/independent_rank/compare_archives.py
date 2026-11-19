"""Check exact overlap against the archived Arb evidence after the new run."""
import csv
import hashlib
import sys
from fractions import Fraction
from pathlib import Path

WITNESSES = ("W1", "W2", "W3")


def read(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def dyadic(mantissa, exponent):
    return Fraction(int(mantissa)) * Fraction(2)**int(exponent)


def old_ball(row):
    midpoint = dyadic(row["mid_man"], row["mid_exp"])
    radius = dyadic(row["rad_man"], row["rad_exp"])
    return midpoint - radius, midpoint + radius


def new_interval(row):
    return tuple(dyadic(row[f"{endpoint}_mantissa"], row[f"{endpoint}_exponent"])
                 for endpoint in ("lower", "upper"))


def overlaps(left, right):
    return max(left[0], right[0]) <= min(left[1], right[1])


def main(stage, observer, own):
    tracked = [stage / "inputs" / name for name in ("cases.tsv", "fixed_settings.tsv", "witnesses.tsv")]
    tracked.append(observer)
    observer_by_witness = {row["witness"]: row for row in read(observer)}
    results = []
    for witness in WITNESSES:
        source = stage / "verification" / ("B04_" + witness)
        tracked += [source / name for name in (
            "candidate_registry.tsv", "branch_registry.tsv", "row_registry.tsv",
            "column_registry.tsv", "primary_128_matrix.tsv", "primary_128_determinant.tsv",
        )]
        original_values = {(int(row["row"]), int(row["column"])): old_ball(row)
                           for row in read(source / "primary_128_matrix.tsv")}
        new_values = {(int(row["row"]), int(row["column"])): new_interval(row)
                      for row in read(own / f"{witness}_matrix.tsv")}
        assert len(new_values) == 361 and len(original_values) == 324
        matches = sum(overlaps(old_value, new_values[key]) for key, old_value in original_values.items())
        sigma_partial_matches = overlaps(old_ball(observer_by_witness[witness]), new_values[(18, 18)])
        exact_zero_column = all(new_values[(i, 18)] == (0, 0) for i in range(18))
        assert matches == 324 and sigma_partial_matches and exact_zero_column
        prior_determinant = old_ball(read(source / "primary_128_determinant.tsv")[0])
        current = read(own / f"{witness}_results.tsv")
        assert len(current) == 2 and {int(item["dimension"]) for item in current} == {18, 19}
        for item in current:
            lower, upper = new_interval(item)
            interval_excludes_zero = lower > 0 or upper < 0
            assert interval_excludes_zero and item["status"] == "FULL_RANK"
            dimension = int(item["dimension"])
            if dimension == 18:
                old_det_matches = overlaps(prior_determinant, (lower, upper))
            else:
                old_sigma = old_ball(observer_by_witness[witness])
                products = [a * b for a in prior_determinant for b in old_sigma]
                old_det_matches = overlaps((min(products), max(products)), (lower, upper))
            assert old_det_matches
            results.append(dict(witness=witness, dimension=dimension, rank=dimension,
                                determinant_sign="positive" if lower > 0 else "negative",
                                residual_norm_upper_rounded=item["residual_norm_upper"],
                                older_matrix_entries_overlap=f"{matches}/324",
                                older_sigma_partial_overlaps=sigma_partial_matches,
                                older_determinant_overlaps=old_det_matches))
    write(own / "summary.tsv", results)
    write(own / "source_sha256.tsv", [
        dict(path=str(path.relative_to(stage.parent)), sha256=hashlib.sha256(path.read_bytes()).hexdigest())
        for path in tracked
    ])


if __name__ == "__main__":
    main(*map(Path, sys.argv[1:]))
