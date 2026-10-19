#!/usr/bin/env python3
"""Build the two-plate synthetic fixture the `test` profile runs on.

The repository carries no real specimen data, so CI exercises the workflow
against plates generated here. They reproduce the structure of a QuantStudio
export exactly: the same sheet name, the same preamble above the header, the
same column set, the Layout A control positions, and `Task == "NTC"` covering
both the no-template controls and the extraction blanks.

Detection is generated with Pf driving it and Pm and Po trailing at a spread
of Cq delays, so the guardband sensitivity arm has something to move at every
value of the grid, as it does in the real data.

Run from the repository root, after the sample sheet:

    python3 tests/fixtures/make_sample_sheet.py
    python3 tests/fixtures/make_plates.py
"""
from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "bin"))
from randomise import slot_wells  # noqa: E402
rng = np.random.default_rng(7)
REPLICATE_CODE='TS005'
TARGETS=('Pf','Pm','Po')
BLOCK={'Pf':0,'Pm':4,'Po':8}
ROWS='ABCDEFGH'
PRE=[['Block Type','96-Well 0.1-mL Block'],['Experiment File Name','synthetic'],
     ['Instrument Type','QuantStudio(TM) 5 System'],['Passive Reference','ROX']]
COLS=['Well','Well Position','Omit','Sample Name','Target Name','Task','Reporter',
      'Quencher','CT','Ct Mean','Ct SD','Quantity','Quantity Mean','Quantity SD',
      'Y-Intercept','R(superscript 2)','Slope','Efficiency',
      'Automatic Ct Threshold','Ct Threshold','Automatic Baseline',
      'Baseline Start','Baseline End','Amp Status','Comments','Cq Conf',
      'HIGHSD','DRNMIN','THOLDFAIL','CQCONF','EXPFAIL']
REP={'Pf':'FAM','Pm':'VIC','Po':'NED'}

def idx(well):
    r,c=well[0],int(well[1:]); return ROWS.index(r)*12+c

def plate(pid, codes, seed):
    rr=np.random.default_rng(seed); rows=[]
    def add(well,name,tgt,task,ct,amp,conf,qty=None):
        rows.append({**{k:None for k in COLS},'Well':idx(well),'Well Position':well,
            'Omit':False,'Sample Name':name,'Target Name':tgt,'Task':task,
            'Reporter':REP[tgt],'Quencher':'NFQ-MGB',
            'CT':'Undetermined' if ct is None else round(float(ct),6),
            'Quantity':qty,'Y-Intercept':40.1,'R(superscript 2)':0.996,
            'Slope':-3.32,'Efficiency':100.1,'Automatic Ct Threshold':True,
            'Automatic Baseline':True,'Amp Status':amp,'Cq Conf':round(conf,6)})
    for t in TARGETS:
        o=BLOCK[t]
        for i in range(1,7):   # STD1-6 duplicate, cols 1-2
            ct=12.4+3.3*(i-1)
            for c in (1,2):
                add(f'{ROWS[i-1]}{o+c}',f'{t}_STD_{i}',t,'STANDARD',
                    ct+rr.normal(0,.08),'Amp',0.99,10.0**(1-i))
        for r_ in ('A','B'):   # NTC
            add(f'{r_}{o+3}',f'NTC_{t}',t,'NTC',None,'No Amp',0.0)
        for r_ in ('C','D'):   # extraction blank
            add(f'{r_}{o+3}',f'ExB_{t}',t,'NTC',None,'No Amp',0.0)
        for r_ in ('G','H'):   # positive control
            add(f'{r_}{o+1}',f'Pos_{t}',t,'Positive control',
                {'Pf':20.3,'Pm':26.9,'Po':27.3}[t]+rr.normal(0,.1),'Amp',0.99)
        for r_ in ('G','H'):   # uninfected host negative control
            add(f'{r_}{o+2}','Neg_Ctl',t,'NEGATIVE',None,'No Amp',0.0)
    for slot,code in enumerate(codes,start=1):
        wells=slot_wells(slot)
        # Pf drives detection; Pm and Po co-detections trail it by a spread
        # of delays, so the guardband sensitivity arm has something to move
        # at every value of the grid -- the structure the real data shows.
        pf_pos = rr.random()<0.65
        pf_cq  = rr.uniform(24,32) if pf_pos else None
        # The inter-plate replicate specimen is forced positive on both
        # plates, with plate-specific jitter, so the concordance stage has a
        # genuinely independent pair to compare.
        if code==REPLICATE_CODE:
            pf_pos, pf_cq = True, 28.0+rr.normal(0,.25)
        for t,well in zip(TARGETS,wells):
            if t=='Pf' and pf_pos:
                add(well,code,t,'UNKNOWN',pf_cq,'Amp',rr.uniform(.9,.99))
            elif t=='Po' and pf_pos and rr.random()<0.75:
                add(well,code,t,'UNKNOWN',pf_cq+rr.uniform(0.2,6.0),'Amp',
                    rr.uniform(.6,.95))
            elif t=='Pm' and pf_pos and rr.random()<0.55:
                add(well,code,t,'UNKNOWN',pf_cq+rr.uniform(0.2,6.0),'Amp',
                    rr.uniform(.6,.95))
            elif t=='Pm' and not pf_pos and rr.random()<0.10:
                add(well,code,t,'UNKNOWN',rr.uniform(30,36),'Amp',
                    rr.uniform(.6,.95))
            else:
                add(well,code,t,'UNKNOWN',None,'No Amp',0.0)
    return pd.DataFrame(rows).sort_values('Well')

def write(path,df):
    from openpyxl import Workbook
    wb=Workbook(); ws=wb.active; ws.title='Results'
    for r in PRE: ws.append(r)
    ws.append([])
    ws.append(COLS)
    for _,row in df[COLS].iterrows():
        ws.append([None if pd.isna(v) else v for v in row.tolist()])
    wb.save(path)

A=[f'TS{i:03d}' for i in range(1,13)]
B=[f'TS{i:03d}' for i in range(13,24)]+['TS005']   # TS005 = inter-plate replicate
OUT = ROOT / 'tests' / 'data' / 'plates'
OUT.mkdir(parents=True, exist_ok=True)
write(OUT / 'Plas_18S_qPCR_R01.xlsx', plate('R01', A, 11))
write(OUT / 'Plas_18S_qPCR_R02.xlsx', plate('R02', B, 12))

print(f'wrote {OUT}/Plas_18S_qPCR_R01.xlsx and R02.xlsx')
