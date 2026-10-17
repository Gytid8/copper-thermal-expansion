"""Small scientific consistency checks; does not create simulated results."""
import numpy as np
from run_study import manifest,source
from analyze_study import reference

jobs=manifest()
assert len(jobs)==24
assert len({j['id'] for j in jobs})==24
prod=[j for j in jobs if j['group']=='production']
assert len(prod)==18 and len({j['seed'] for j in prod})==18
for t in [300,400,500,600,700,800]:
    assert sum(j['temperature']==t for j in prod)==3
for j in jobs:
    assert round(.1/j['dt'])*j['dt']==.1
    assert round(j['sample_ps']/j['dt'])*j['dt']==j['sample_ps']
    s=source(j)
    assert 'npt temp' in s and 'iso 0 0 1' in s
    assert 'format " %.12g"' in s
ref,alpha=reference()
assert ref.loc[ref.temperature_K==300,'strain_from_300'].iloc[0]==0
assert np.isclose(alpha,(9420-117)*1e-6/(1+117e-6)/500,rtol=1e-12)
assert np.isclose(alpha*1e6,18.6038233527,rtol=1e-8)
assert np.all(np.diff(ref.strain_from_300)>0)
print('PASS: 24 jobs, unique seeds, physical durations, spaced output, 293→300 K rebasing')
