import numpy as np
import glob, re, sys

def load_chain(prefix):
    files = sorted(f for f in glob.glob(prefix + '.*.txt') if re.search(r'\.\d+\.txt$', f))
    if not files:
        return None, None, None
    blocks = []
    for f in files:
        with open(f) as fh:
            header = fh.readline()
        cols = header.split()[1:]  # 跳过 '#'
        data = np.loadtxt(f, comments='#')
        blocks.append(data)
    d = np.concatenate(blocks)
    wi, wai = cols.index('w'), cols.index('wa')
    return d, wi, wai

def analyze(prefix):
    d, wi, wai = load_chain(prefix)
    if d is None:
        return None
    w, wa = d[:, wi], d[:, wai]
    cov = np.cov(np.vstack([w, wa]))
    delta = np.array([w.mean() + 1, wa.mean()])
    chi2 = delta @ np.linalg.inv(cov) @ delta
    sig = np.sqrt(chi2)
    return len(d), w.mean(), w.std(), wa.mean(), wa.std(), \
           np.corrcoef(w, wa)[0, 1], sig

if __name__ == '__main__':
    dataset = sys.argv[1] if len(sys.argv) > 1 else 'dovekie'
    base = rf'D:\research\cosmology\truncation_sn\{dataset}_full_planck\fp_{dataset}'
    zcuts = ['z0p00', 'z0p02', 'z0p05', 'z0p06', 'z0p07', 'z0p08', 'z0p09',
             'z0p10', 'z0p11', 'z0p12', 'z0p14', 'z0p15', 'z0p20']
    print(f"{'zcut':<6} {'N':>8} {'w0':>10} {'+-':>6} {'wa':>10} {'+-':>6} {'corr':>7} {'sig2D':>7}")
    for z in zcuts:
        r = analyze(base + '_' + z)
        if r is None:
            print(f"{z:<6}  missing")
            continue
        N, w0, e0, wa, ea, corr, sig = r
        print(f"{z:<6} {N:>8} {w0:>10.4f} {e0:>6.4f} {wa:>10.4f} {ea:>6.4f} {corr:>7.3f} {sig:>7.2f}")
