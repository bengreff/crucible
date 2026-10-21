#!/usr/bin/env python3
"""Check an engine-spec file (engine_spec.schema.json) and screen it for outcome values.

  check_specs.py FILE                 schema and consistency check, then the usability table (inputs only)
  check_specs.py FILE --counts        the outcome screen, count-only: hits per keyword, no text printed
  check_specs.py FILE --masked        every hit's path and its text with each digit run shown as '#'
  check_specs.py FILE --redact PATH   replace the string at PATH unread (a num field's value: null, missing)

The screen reads every free-text string (notes, original, derivation, locator, citations, labels, names,
flag reasons, URLs). Enum values are not screened. Stdlib only: the validator covers the schema subset the
schema uses ($ref, allOf, type, enum, const, required, properties, additionalProperties, items, pattern).
"""
import json, re, sys
from pathlib import Path

SCHEMA = json.loads((Path(__file__).parent / "engine_spec.schema.json").read_text())

OUTCOME = {
    "isp": r"\bisp\b|specific impulse|\bi_?sp\b",
    "thrust": r"thrust(?![ _-]?(optimi[sz]ed|chamber|vector))",
    "c_star": r"c\s?\*|c-star|\bcstar\b|characteristic velocity",
    "exhaust_velocity": r"exhaust velocity|\bv_?e\b|\bc_?eff\b",
    "thrust_coefficient": r"thrust coefficient|\bc_?f\b",
    "force_units": r"\b(k?n|mn|k?lbf|klb|kgf|tf|tonnes?|newtons?)\b",
    "mass_flow": r"mass flow|flow ?rate|\bm_?dot\b|kg/s|lbm?/s|lb/sec|\bkg/sec\b",
    "efficiency": r"efficien",
    "seconds": r"\bseconds?\b|\bsec\b(?!\.)|\d\s*s\b",
    "velocity_units": r"\bm/s\b|\bft/s\b|\bkm/s\b",
    "exit_pressure": r"exit pressure|\bp_?e\b",
}
FREE_TEXT = {"note", "original", "derivation", "locator", "citation", "label", "name", "notes", "reason", "url",
             "maker", "country", "quantity_label"}
ENUM_VALUE_PATHS = ("cycle.type", "cycle.turbine_exhaust", "operating_point.pressure_station",
                    "operating_point.mixture_ratio_basis", "geometry.contour_type", "geometry.chamber_cooling",
                    "injector.type")


def resolve(ref):
    node = SCHEMA
    for part in ref.lstrip("#/").split("/"):
        node = node[part]
    return node


def validate(x, s, path, errs):
    if "$ref" in s:
        validate(x, resolve(s["$ref"]), path, errs)
    for sub in s.get("allOf", []):
        validate(x, sub, path, errs)
    if "const" in s and x != s["const"]:
        errs.append(f"{path}: must be {s['const']!r}")
    if "enum" in s and x not in s["enum"]:
        errs.append(f"{path}: {x!r} not in {s['enum']}")
    t = s.get("type")
    if t:
        ts = t if isinstance(t, list) else [t]
        ok = {"object": dict, "array": list, "string": str, "null": type(None)}
        good = any((tt == "number" and isinstance(x, (int, float)) and not isinstance(x, bool))
                   or (tt in ok and isinstance(x, ok[tt])) for tt in ts)
        if not good:
            errs.append(f"{path}: type {type(x).__name__}, wants {t}")
            return
    if isinstance(x, str) and "pattern" in s and not re.search(s["pattern"], x):
        errs.append(f"{path}: {x!r} fails {s['pattern']}")
    if isinstance(x, dict):
        props = s.get("properties", {})
        for r in s.get("required", []):
            if r not in x:
                errs.append(f"{path}: missing '{r}'")
        if s.get("additionalProperties") is False:
            for k in x:
                if k not in props:
                    errs.append(f"{path}: key '{k}' not allowed")
        for k, sub in props.items():
            if k in x and sub:
                validate(x[k], sub, f"{path}.{k}", errs)
    if isinstance(x, list) and "items" in s:
        for i, v in enumerate(x):
            validate(v, s["items"], f"{path}[{i}]", errs)


def fields(e):
    """(dotted path, field object) for every value/provenance field of an engine."""
    for sect in ("propellants", "cycle", "operating_point", "geometry", "injector", "temperatures"):
        for k, f in e.get(sect, {}).items():
            if isinstance(f, dict) and "provenance" in f:
                yield f"{sect}.{k}", f


def consistency(e, i, errs):
    ids = {s["id"] for s in e.get("sources", [])}
    for p, f in fields(e):
        pv, here = f.get("provenance"), f"engines[{i}]({e.get('id')}).{p}"
        if pv in ("printed", "digitised") and not (f.get("source") and f.get("locator")):
            errs.append(f"{here}: {pv} needs source and locator")
        if pv == "derived" and not f.get("derivation"):
            errs.append(f"{here}: derived needs derivation")
        if pv == "missing" and f.get("value") is not None:
            errs.append(f"{here}: missing must have value null")
        if pv != "missing" and f.get("value") is None:
            errs.append(f"{here}: value null but provenance {pv}")
        if f.get("source") and f["source"] not in ids:
            errs.append(f"{here}: source '{f['source']}' not in sources")
    for j, o in enumerate(e.get("outcome_pointers", [])):
        if o["source"] not in ids:
            errs.append(f"engines[{i}].outcome_pointers[{j}]: source '{o['source']}' not in sources")


def strings(x, path=""):
    if isinstance(x, dict):
        for k, v in x.items():
            p = f"{path}.{k}" if path else k
            if isinstance(v, str) and (k in FREE_TEXT or (k == "value" and not path.endswith(ENUM_VALUE_PATHS))):
                yield p, v
            else:
                yield from strings(v, p)
    elif isinstance(x, list):
        for j, v in enumerate(x):
            yield from strings(v, f"{path}[{j}]")


def hits(doc):
    for p, s in strings(doc):
        kws = [k for k, r in OUTCOME.items() if re.search(r, s, re.I)]
        if kws:
            yield p, s, kws


def have(e, p):
    sect, k = p.split(".")
    f = e.get(sect, {}).get(k)
    return bool(f) and f.get("provenance") != "missing" and f.get("value") not in (None, "unknown", "unspecified")


TIER0_ISP = ["propellants.oxidizer", "propellants.fuel", "operating_point.chamber_pressure_mpa",
             "operating_point.mixture_ratio", "geometry.expansion_ratio"]
TIER0_THRUST = TIER0_ISP + ["geometry.throat_diameter_m"]
TIER1 = TIER0_THRUST + ["geometry.contraction_ratio", "geometry.contour_type", "injector.type",
                        "temperatures.fuel_injection_k", "temperatures.oxidizer_injection_k"]


def usability(e):
    out = {}
    for name, req in (("tier0_isp", TIER0_ISP), ("tier0_thrust", TIER0_THRUST), ("tier1", TIER1)):
        miss = [p.split(".")[1] for p in req if not have(e, p)]
        if name == "tier1" and not (have(e, "geometry.chamber_length_m") or have(e, "geometry.characteristic_length_m")):
            miss.append("chamber_length|L*")
        out[name] = miss
    return out


def walk_to(doc, path):
    keys = re.findall(r"[^.\[\]]+|\[\d+\]", path)
    node = doc
    for k in keys[:-1]:
        node = node[int(k[1:-1])] if k.startswith("[") else node[k]
    last = keys[-1]
    return node, (int(last[1:-1]) if last.startswith("[") else last)


def main():
    a = sys.argv[1:]
    if not a:
        sys.exit(__doc__)
    fn = Path(a[0])
    doc = json.loads(fn.read_text())
    if "--counts" in a:
        n = {k: 0 for k in OUTCOME}
        for _, _, kws in hits(doc):
            for k in kws:
                n[k] += 1
        print(json.dumps(n))
        return
    if "--masked" in a:
        for p, s, kws in hits(doc):
            print(f"{p} [{','.join(kws)}]: {re.sub(r'[0-9][0-9.,]*', '#', s)}")
        return
    if "--redact" in a:
        node, k = walk_to(doc, a[a.index("--redact") + 1])
        if k == "value" and isinstance(node, dict) and "provenance" in node:
            node.update(value=None, provenance="missing", note="removed unread by the outcome screen",
                        original=None, derivation=None, locator=None)
        else:
            node[k] = "[removed unread by the outcome screen]"
        fn.write_text(json.dumps(doc, indent=1, ensure_ascii=False) + "\n")
        print("redacted", a[a.index("--redact") + 1])
        return
    errs = []
    validate(doc, SCHEMA, "$", errs)
    for i, e in enumerate(doc.get("engines", [])):
        consistency(e, i, errs)
    ids = [e.get("id") for e in doc.get("engines", [])]
    errs += [f"duplicate id {x}" for x in sorted({x for x in ids if ids.count(x) > 1})]
    print(f"{len(ids)} engines, {len(errs)} errors")
    for x in errs:
        print("  " + x)
    print(f"{'id':32s} {'family':13s} tier0_isp tier0_thrust tier1  missing for tier1")
    for e in doc.get("engines", []):
        u = usability(e)
        mark = lambda m: "yes" if not m else "no"
        print(f"{e.get('id', '?'):32s} {e.get('family', '?'):13s} {mark(u['tier0_isp']):9s} "
              f"{mark(u['tier0_thrust']):12s} {mark(u['tier1']):5s}  {','.join(u['tier1'])}")
    sys.exit(1 if errs else 0)


if __name__ == "__main__":
    main()
