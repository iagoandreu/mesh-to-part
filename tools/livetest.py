"""Live self-test: drives the running add-in through the dev bridge and prints
everything. Waits for the bridge to come up first.

    py -3.10 livetest.py [wait_seconds]
"""

import json
import sys
import time

import mtp


def show(title, res):
    print('\n' + '=' * 70)
    print('## ' + title)
    print('=' * 70)
    if not res.get('ok'):
        print('FAILED: ' + str(res.get('error'))[:4000])
        return None
    r = res.get('result', {})
    text = r.pop('text', None) if isinstance(r, dict) else None
    print(json.dumps(r, indent=2, default=str)[:4000])
    if text:
        print('\n' + text)
    return r


def wait_for_bridge(limit):
    t0 = time.time()
    while time.time() - t0 < limit:
        if mtp.alive():
            print('bridge is up after %.0f s' % (time.time() - t0))
            return True
        time.sleep(2.0)
    return False


def main():
    limit = float(sys.argv[1]) if len(sys.argv) > 1 else 900.0
    if not wait_for_bridge(limit):
        print('bridge never came up - the add-in has not been re-run')
        return 1

    show('ping - what is open', mtp.send('ping'))
    show('panels - where the button went', mtp.send('panels'))
    show('probe - which conversions this licence allows', mtp.send('probe'))
    show('analyse - gaps in the real STL', mtp.send('analyse'))

    # --- the real STL, cloned into a throw-away document -------------------- #
    cl = show('clone the STL into a scratch document', mtp.send('clone_to_new_doc'))
    if cl and not cl.get('error'):
        show('analyse the clone', mtp.send('analyse'))
        for name, opts in (
            ('PRISMATIC on the real STL',
             {'method_key': 'prismatic', 'max_error_mm': 0.1}),
            ('FACETED on the real STL',
             {'method_key': 'faceted', 'max_error_mm': 0.1, 'escalate': False}),
            ('ORGANIC on the real STL (expected to be refused cleanly)',
             {'method_key': 'organic', 'max_error_mm': 0.1}),
        ):
            show(name, mtp.send('convert', {'options': opts}, timeout=900))
            show('reset scratch doc', mtp.send('close_doc'))
            show('re-clone', mtp.send('clone_to_new_doc'))
        show('capture', mtp.send('capture'))
        show('close scratch doc', mtp.send('close_doc'))

    # --- synthetic sphere with a known hole --------------------------------- #
    show('make synthetic test doc', mtp.send('make_test_doc'))
    show('analyse the sphere (hole should be ~7.765 mm)', mtp.send('analyse'))
    show('PRISMATIC on the sphere, gap limit 8 mm',
         mtp.send('convert', {'options': {'method_key': 'prismatic',
                                          'max_error_mm': 0.1, 'max_gap_mm': 8.0}}))
    show('capture sphere', mtp.send('capture'))
    show('close sphere doc', mtp.send('close_doc'))

    show('ping - back where we started', mtp.send('ping'))
    print('\nLIVE TEST COMPLETE')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
