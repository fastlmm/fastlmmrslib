"""Differential qualification against installed legacy artifacts and saved inputs."""
import json
import subprocess
from pathlib import Path

import numpy as np
from fastlmmclib.quadform import qf

root = Path(__file__).resolve().parent
evidence = Path('/tmp/fastlmm-maintenance-review')
cases = []
for line in (evidence / 'gx2-evaluation.jsonl').read_text().splitlines():
    row = json.loads(line)
    if 'gx2' in row:
        cases.append(dict(name=row['name'], x=row['x'], w=row['w'], k=row['k'],
                          nc=row['nc'], sigma=row['sigma'], acc=1e-10, lim=1_000_000))
for index, line in enumerate((evidence / 'gx2-real-inputs.jsonl').read_text().splitlines()):
    row = json.loads(line)
    cases.append(dict(name=f'fastlmm_{index}', x=row['cutoff'], w=row['weights'],
                      k=[1]*len(row['weights']), nc=[0.0]*len(row['weights']),
                      sigma=0.0, acc=1e-7, lim=1_000_000))
rng = np.random.default_rng(20261113)
for index in range(100):
    n = int(rng.choice([1, 2, 10, 50]))
    w = rng.normal(size=n)
    k = rng.integers(0, 5, size=n).astype('int32')
    nc = rng.uniform(0.0, 3.0, size=n)
    sigma = float(rng.choice([0.0, 0.5, 2.0]))
    mean = float(np.sum(w*(k+nc)))
    sd = float(np.sqrt(np.sum(w*w*(2*k+4*nc))+sigma*sigma))
    cases.append(dict(name=f'generic_{index}', x=mean+float(rng.choice([-3,0,3]))*sd,
                      w=w.tolist(), k=k.tolist(), nc=nc.tolist(), sigma=sigma,
                      acc=1e-7, lim=1_000_000))

lines = []
for c in cases:
    values = [len(c['w']), c['x'], c['sigma'], c['lim'], c['acc'], *c['w'], *c['k'], *c['nc']]
    lines.append(' '.join(str(v) for v in values))
(root/'comparison-cases.json').write_text(json.dumps(cases))
with (root/'comparison-progress.log').open('w') as progress:
    run = subprocess.run([str(root/'target/release/examples/probe')],
                         input='\n'.join(lines)+'\n', text=True,
                         stdout=subprocess.PIPE, stderr=progress, check=True, timeout=30)
answers = run.stdout.splitlines()
assert len(answers) == len(cases)
results = []
for case, answer in zip(cases, answers):
    native, fault, trace = qf(case['x'], np.array(case['w'], dtype='float64'),
        dof=np.array(case['k'], dtype='int32'), noncentrality=np.array(case['nc'], dtype='float64'),
        sigma=case['sigma'], acc=case['acc'], lim=case['lim'])
    fields = answer.split()
    rust_fault = int(fields[0])
    rust = None if fields[1]=='none' else float(fields[1])
    rust_trace = np.array([float(x) for x in fields[2:]])
    assert len(rust_trace)==7
    results.append(dict(name=case['name'], qfc=float(native), fault=int(fault),
                        rust=rust, rust_fault=rust_fault,
                        probability_difference=None if rust is None else rust-float(native),
                        max_trace_difference=float(np.max(np.abs(rust_trace-trace)))))
output = root/'comparison-results.json'
output.write_text(json.dumps(results, indent=2)+'\n')
mismatched = [r for r in results if r['rust_fault']!=r['fault'] or
              (r['fault'] in [0,2] and (r['rust'] is None or abs(r['probability_difference'])>1e-12)) or
              (r['fault'] not in [0,2] and r['rust'] is not None)]
print(json.dumps(dict(cases=len(results), mismatches=len(mismatched),
    max_probability_difference=max(abs(r['probability_difference']) for r in results if r['rust'] is not None),
    max_trace_difference=max(r['max_trace_difference'] for r in results),
    mismatch_details=mismatched), indent=2))
assert not mismatched
