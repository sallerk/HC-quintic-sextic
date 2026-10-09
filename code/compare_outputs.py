"""compare_outputs.py -- compare the JSON files of a rerun with the reference files in outputs/.
Keys named "seconds" (run times) are ignored; floating-point numbers are compared with relative tolerance 1e-9;
everything else must be equal.
usage (from the rerun folder):  python ../code/compare_outputs.py ../outputs .
Exit code 0 if every reference .json file has an equal counterpart, 1 otherwise.
"""
import json, math, os, sys


def diff(a, b, path=''):
    if isinstance(a, float) or isinstance(b, float):
        try:
            return [] if math.isclose(float(a), float(b), rel_tol=1e-9, abs_tol=1e-12) else [path]
        except (TypeError, ValueError):
            return [path]
    if type(a) != type(b):
        return [path]
    if isinstance(a, dict):
        out = [path + '/' + k for k in sorted(set(a) ^ set(b))]
        for k in sorted(set(a) & set(b)):
            if k != 'seconds':
                out += diff(a[k], b[k], path + '/' + k)
        return out
    if isinstance(a, list):
        if len(a) != len(b):
            return [path + ' (length %d vs %d)' % (len(a), len(b))]
        out = []
        for i, (x, y) in enumerate(zip(a, b)):
            out += diff(x, y, path + '[%d]' % i)
        return out
    return [] if a == b else [path]


def main(ref, new):
    ok = True
    for name in sorted(f for f in os.listdir(ref) if f.endswith('.json')):
        p = os.path.join(new, name)
        if not os.path.exists(p):
            print('%-22s MISSING in %s' % (name, new))
            ok = False
            continue
        d = diff(json.load(open(os.path.join(ref, name))), json.load(open(p)))
        print('%-22s %s' % (name, 'equal' if not d else '%d differences, first: %s' % (len(d), d[:5])))
        ok &= not d
    print('ALL EQUAL' if ok else 'DIFFERENCES FOUND')
    return ok


if __name__ == '__main__':
    sys.exit(0 if main(sys.argv[1], sys.argv[2]) else 1)
