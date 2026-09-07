"""Translate the pinned ParaEMT 9-bus data into GridKit's continuous EMT model."""
import csv
import cmath
import json
import math
import copy
from pathlib import Path

import numpy as np
from scipy.optimize import root

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'gridkit'


def table(name):
    with (ROOT / name).open() as stream:
        return {row[0]: row[1:] for row in csv.reader(stream)}


def diagonal(value):
    return [[value if i == j else 0.0 for j in range(3)] for i in range(3)]


def main():
    OUT.mkdir(exist_ok=True)
    pf = json.loads((ROOT / 'pfd_9_1_1.json').read_text())
    ec = json.loads((ROOT / 'results/governor_step/dt25us/machine_parameters.json').read_text())
    gen, exc, gov, pss = (table(name + '.csv') for name in ('gen', 'exc', 'gov', 'pss'))
    omega, voltage, power = pf['ws'], 230e3, pf['basemva'] * 1e6
    zb = voltage**2 / power
    # The published JSON rounds voltages and dispatch independently. Resolve
    # the passive constant-impedance network with its original load impedances,
    # fixed generator voltage magnitudes, and fixed P at generators 2 and 3.
    # This supplies consistent inductor currents without adding capacitance.
    admittance = np.zeros((9, 9), dtype=complex)
    for kind in ('line', 'xfmr'):
        for i, (a, b, raw) in enumerate(zip(pf[kind+'_from'], pf[kind+'_to'], pf[kind+'_RX'])):
            a, b = a-1, b-1
            y = 1/complex(raw)
            admittance[a,a] += y; admittance[b,b] += y
            admittance[a,b] -= y; admittance[b,a] -= y
            if kind == 'line':
                admittance[a,a] += .5j*pf['line_chg'][i]
                admittance[b,b] += .5j*pf['line_chg'][i]
    for i, bus in enumerate(pf['load_bus']):
        admittance[bus-1,bus-1] += complex(pf['load_MW'][i],-pf['load_Mvar'][i])/pf['basemva']/pf['bus_Vm'][bus-1]**2
    transfer = -np.linalg.solve(admittance[3:,3:], admittance[3:,:3])
    reduced = admittance[:3,:3] + admittance[:3,3:] @ transfer
    magnitudes = np.array(pf['bus_Vm'][:3])
    def generation(angles):
        vg = magnitudes*np.exp(1j*np.r_[0.,angles])
        return vg, vg*np.conj(reduced @ vg)
    solved = root(lambda angles: generation(angles)[1].real[1:] - np.array(pf['gen_MW'][1:])/pf['basemva'], pf['bus_Va'][1:3], tol=1e-11)
    if not solved.success and np.max(np.abs(solved.fun)) > 1e-12:
        raise RuntimeError(solved.message)
    vg, sg = generation(solved.x)
    balanced = np.r_[vg, transfer @ vg]
    original = np.array(pf['bus_Vm'])*np.exp(1j*np.array(pf['bus_Va']))
    adjustment = {'method': 'Schur reduction of the original constant-Z network; fixed generator Vm, P2, P3 and angle1',
                  'max_voltage_phasor_adjustment_pu': float(np.max(np.abs(balanced-original))),
                  'bus_Vm_pu': abs(balanced).tolist(), 'bus_Va_rad': np.angle(balanced).tolist(),
                  'gen_MW': (sg.real*pf['basemva']).tolist(), 'gen_Mvar': (sg.imag*pf['basemva']).tolist(),
                  'gen_MW_adjustment': (sg.real*pf['basemva']-pf['gen_MW']).tolist(),
                  'gen_Mvar_adjustment': (sg.imag*pf['basemva']-pf['gen_Mvar']).tolist(),
                  'max_passive_bus_current_mismatch_pu': float(np.max(np.abs((admittance @ balanced)[3:]))) }
    (OUT/'initialization.json').write_text(json.dumps(adjustment,indent=2)+'\n')
    v = {bus: math.sqrt(2/3)*voltage*balanced[bus-1] for bus in pf['bus_num']}
    rotation = [cmath.rect(1, shift) for shift in (0, -2*math.pi/3, 2*math.pi/3)]
    devices, signals, seed, reduction = [], [], {}, {}
    state = {'header': {'version': 1, 'time': 0.0}, 'buses': {}, 'devices': {}}

    def initial(name, phasors):
        seed[name] = {'y': [z.real for z in phasors], 'yp': [(1j*omega*z).real for z in phasors]}

    for bus in pf['bus_num']:
        name = f'bus_{bus}'
        phases = [v[bus] * r for r in rotation]
        initial(name, phases)
        state['buses'][name] = dict(zip(('va', 'vb', 'vc'), seed[name]['y']))
        devices.append({'id': name, 'class': 'Bus', 'mon': ['va', 'vb', 'vc']})

    for i, bus in enumerate(pf['gen_bus']):
        signals.extend({'id': f'{name}_{bus}'} for name in ('speed', 'pm', 'efd', 'pss'))
        signals.append({'id': f'pref_{bus}'})
        params = {'N': 3, 'S': pf['gen_MVA_base'][i]*1e6, 'V': voltage, 'f': omega/(2*math.pi),
                  'H': float(gen['H'][i]), 'F': float(gen['D'][i]), 'S10': float(gen['S(1.0)'][i]), 'S12': float(gen['S(1.2)'][i])}
        params.update({key: ec['ec_'+val][i] for key, val in [('Ll','Ll'),('Lmd','Lad'),('Lmq','Laq'),('L0','L0'),('Rs','Ra')]})
        params.update({key: ec['ec_'+total][i]-ec['ec_'+mutual][i] for key,total,mutual in [('Llfd','Lffd','Lad'),('Ll1d','L11d','Lad'),('Ll1q','L11q','Laq'),('Ll2q','L22q','Laq')]})
        params.update({key: ec['ec_'+key][i]/omega for key in ('Rfd','R1d','R1q','R2q')})
        ld, lq = params['Lmd'], params['Lmq']
        md = np.array([[-ld-params['Ll'],ld,ld],[-ld,ld+params['Llfd'],ld],[-ld,ld,ld+params['Ll1d']]])
        mq = np.array([[-lq-params['Ll'],lq,lq],[-lq,lq+params['Ll1q'],lq],[-lq,lq,lq+params['Ll2q']]])
        reduction[str(bus)] = {'d': np.linalg.inv(md)[0].tolist(), 'q': np.linalg.inv(mq)[0].tolist(),
                               'L0': params['L0'], 'omega': omega}
        assert params['S10'] == params['S12'] == params['F'] == 0
        devices.append({'id': f'gen_{bus}', 'class': 'Machine', 'params': params,
                        'inputs': {'bus': f'bus_{bus}', 'pm': f'pm_{bus}', 'efd': f'efd_{bus}'},
                        'outputs': {'speed': f'speed_{bus}'}, 'mon': ['omega','efd','p','q']})
        state['devices'][f'gen_{bus}'] = {'p': float(sg[i].real*power), 'q': float(sg[i].imag*power)}
        # The ParaEMT bus measurement sheet has te=0.02 s on every bus.
        ep = {'V': voltage, 'Tr': .02, 'Ta': float(exc['TA_o_TB'][i])*float(exc['TB'][i])}
        ep.update({key: float(exc[src][i]) for key,src in [('Tb','TB'),('Te','TE'),('K','K'),('Efdmin','Emin'),('Efdmax','Emax')]})
        devices.append({'id': f'exc_{bus}', 'class': 'SEXS-PTI', 'params': ep,
                        'inputs': {'bus': f'bus_{bus}', 'vs': f'pss_{bus}'}, 'outputs': {'efd': f'efd_{bus}'}, 'mon': ['vts']})
        gp = {'S': params['S'], 'Trate': params['S']/1e6}
        gp.update({key: float(gov[src][i]) for key,src in [('R','R'),('T1','T1'),('T2','T2'),('T3','T3'),('At','Ambient temperature load limit'),('Kt','KT'),('Vmax','VMAX'),('Vmin','VMIN'),('Dturb','Dturb')]})
        devices.append({'id': f'gov_{bus}', 'class': 'GASTPTI', 'params': gp, 'inputs': {'speed': f'speed_{bus}', 'pref': f'pref_{bus}'}, 'outputs': {'pmech': f'pm_{bus}'}, 'mon': ['pmech']})
        pp = {key: float(pss[key][i]) for key in ('A1','A2','A3','A4','A5','A6','T1','T2','T3','T4','T5','T6')}
        pp.update({key: float(pss[src][i]) for key,src in [('Ks','KS'),('Lsmin','LSMIN'),('Lsmax','LSMAX'),('Vcl','VCL'),('Vcu','VCU')]})
        devices.append({'id': f'stabilizer_{bus}', 'class': 'IEEEST', 'params': pp, 'inputs': {'speed': f'speed_{bus}'}, 'outputs': {'output': f'pss_{bus}'}, 'mon': ['vss']})

    for kind in ('line', 'xfmr'):
        for i, (a,b,raw) in enumerate(zip(pf[kind+'_from'], pf[kind+'_to'], pf[kind+'_RX'])):
            if kind == 'xfmr': assert pf['xfmr_k'][i] == 1
            z = complex(raw)*zb
            c = pf['line_chg'][i]/(omega*zb) if kind == 'line' else 0.0
            name = f'{kind}_{a}_{b}'
            devices.append({'id': name, 'class': 'LineLumped', 'params': {'N': 3, 'K': 3, 'conductors': [1,2,3], 'dx': 1.0, 'Rp': diagonal(z.real), 'Lp': diagonal(z.imag/omega), 'Gp': diagonal(0.0), 'Cp': diagonal(c)}, 'inputs': {'bus1': f'bus_{a}', 'bus2': f'bus_{b}'}})
            initial(name, [(v[a]-v[b])/z*r for r in rotation] + [-1j*omega*c/2*v[a]*r for r in rotation] + [-1j*omega*c/2*v[b]*r for r in rotation])
    for i,bus in enumerate(pf['load_bus']):
        index = pf['bus_num'].index(bus)
        z = pf['bus_Vm'][index]**2 / complex(pf['load_MW'][i], -pf['load_Mvar'][i]) * pf['basemva'] * zb
        assert z.imag > 0
        name = f'load_{bus}'
        devices.append({'id': name, 'class': 'LoadZ', 'params': {'N': 3, 'R': diagonal(z.real), 'L': diagonal(z.imag/omega)}, 'inputs': {'bus': f'bus_{bus}'}})
        initial(name, [-v[bus]/z*r for r in rotation])

    case = {'header': {'case_name': 'ParaEMT 9-bus', 'case_description': 'Three GENROU equivalent circuits with SEXS, GAST and IEEEST', 'case_comments': 'All voltages referred to 230 kV; original per-unit power bases retained.'}, 'devices': devices, 'signals': signals}
    for name, content in [('9bus.case.json',case),('9bus.state.json',state),('network.state.json',seed),('index_reduction.json',reduction)]:
        (OUT/name).write_text(json.dumps(content,indent=2)+'\n')
    for label,tol in [('tol1e-7',1e-7),('tol1e-8',1e-8),('tol1e-9',1e-9)]:
        study = {'system_model_file': '9bus.case.json', 'state_file': '9bus.state.json', 'dt_monitor': .0005, 'tmax': 3.0, 'rel_tol': tol, 'abs_tol': tol, 'max_steps': 1000000, 'mu': 50000.0, 'consistent_ic_type': 'ya_ydp', 'output_file': 'mon.csv', 'state_output_file': 'state.csv'}
        (OUT/(label+'.solver.json')).write_text(json.dumps(study,indent=2)+'\n')
        study.update(system_model_file='9bus-trip.case.json',state_file='9bus-trip.state.json')
        (OUT/('trip-'+label+'.solver.json')).write_text(json.dumps(study,indent=2)+'\n')
    trip_case,trip_state=copy.deepcopy(case),copy.deepcopy(state)
    terminal='machine_bus_1'
    for device in trip_case['devices']:
        if device['id'] in ('gen_1','exc_1'): device['inputs']['bus']=terminal
    trip_case['devices'].append({'id':terminal,'class':'Bus','mon':['va','vb','vc']})
    trip_state['buses'][terminal]=copy.deepcopy(state['buses']['bus_1'])
    (OUT/'9bus-trip.case.json').write_text(json.dumps(trip_case,indent=2)+'\n')
    (OUT/'9bus-trip.state.json').write_text(json.dumps(trip_state,indent=2)+'\n')


if __name__ == '__main__':
    main()
