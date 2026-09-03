"""Round 2: verify the fixes the first live run exposed.

    graceful failure of a method the geometry cannot support
    clean refusal of a method the licence does not have
    the faceted cap
    no leftover [work] bodies afterwards
"""

import json

import mtp


def show(title, res, keep=('warnings', 'attempts', 'ladder', 'rung_used', 'fell_back',
                           'methods_available', 'faceted_capped', 'result',
                           'faces_out', 'solid_out', 'work_body', 'secs',
                           'reduction_dropped', 'tri_in', 'tri_reduced')):
    print('\n' + '=' * 70)
    print('## ' + title)
    print('=' * 70)
    if not res.get('ok'):
        print('ERROR RETURNED:')
        print(str(res.get('error'))[:1500])
        return None
    r = res.get('result', {})
    if isinstance(r, dict):
        r.pop('text', None)
        if 'deviation' in r and r['deviation']:
            d = r['deviation']
            r['dev_max_mm'] = None if not d.get('out') else round(d['out']['max'], 5)
            r['dev_back_mm'] = None if not d.get('back') else round(d['back']['max'], 5)
        slim = {k: v for k, v in r.items() if k in keep or k.startswith('dev_')}
        print(json.dumps(slim or r, indent=2, default=str)[:3000])
    else:
        print(r)
    return r


def clone():
    return mtp.send('clone_to_new_doc')


def main():
    if not mtp.alive():
        print('bridge is down')
        return 1

    show('start clean', mtp.send('purge_work'))

    clone()
    show('PRISMATIC on the calcaneus - must fall back, not crash',
         mtp.send('convert', {'options': {'method_key': 'prismatic',
                                          'max_error_mm': 0.5}}))
    show('leftovers after fallback', mtp.send('purge_work'))
    mtp.send('close_doc')

    clone()
    show('ORGANIC with fallback allowed - must convert as Faceted',
         mtp.send('convert', {'options': {'method_key': 'organic',
                                          'max_error_mm': 0.5}}))
    mtp.send('close_doc')

    clone()
    show('ORGANIC with fallback refused - must fail with one clear line',
         mtp.send('convert', {'options': {'method_key': 'organic',
                                          'fallback_faceted': False}}))
    show('leftovers after a clean refusal (must be empty)', mtp.send('purge_work'))
    mtp.send('close_doc')

    clone()
    show('FACETED with a 3000-face cap on 10876 triangles',
         mtp.send('convert', {'options': {'method_key': 'faceted',
                                          'faceted_cap': 3000,
                                          'max_error_mm': 2.0,
                                          'escalate': False}}))
    show('capture', mtp.send('capture', {'path': None}))
    mtp.send('close_doc')

    show('synthetic sphere doc', mtp.send('make_test_doc'))
    show('PRISMATIC on the sphere - regression check',
         mtp.send('convert', {'options': {'method_key': 'prismatic',
                                          'max_error_mm': 0.1,
                                          'max_gap_mm': 8.0}}))
    mtp.send('close_doc')

    show('final state', mtp.send('ping'))
    print('\nROUND 2 COMPLETE')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
