import json
import sys

from dimorphite_dl import protonate_smiles


def main():
    payload = json.load(sys.stdin)
    ph = float(payload["ph"])
    smiles_list = list(payload["smiles"])

    results = {}
    for smiles in smiles_list:
        try:
            results[smiles] = protonate_smiles(
                smiles,
                ph_min=ph,
                ph_max=ph,
                precision=0.0,
                max_variants=32,
            )
        except Exception:
            results[smiles] = []

    json.dump(results, sys.stdout)


if __name__ == "__main__":
    main()
