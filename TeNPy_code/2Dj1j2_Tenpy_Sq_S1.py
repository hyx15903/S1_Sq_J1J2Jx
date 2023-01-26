import numpy as np
import tenpy.linalg.np_conserved as npc
import pickle
import gzip

from tenpy.networks.mps import MPS
from tenpy.models.spins_nnn import SpinChainNNN2
from tenpy.algorithms import dmrg

from tenpy.models.lattice import Honeycomb
from tenpy.models.model import CouplingMPOModel
from tenpy.models.model import MultiCouplingModel
from tenpy.tools.params import get_parameter
from tenpy.networks.site import SpinSite



class Chiral_term_for_Square(CouplingMPOModel, MultiCouplingModel):
    def __init__(self, model_params):
        CouplingMPOModel.__init__(self, model_params)

    def init_sites(self, model_params):
        S = model_params.get('S', 1)
        conserve = model_params.get('conserve', 'Sz')
        site = SpinSite(S, conserve)
        return site

    def init_terms(self, model_params):
        Jx = model_params.get('Jx', 1.)
        Jy = model_params.get('Jy', 1.)
        Jz = model_params.get('Jz', 1.)
        Jxp = model_params.get('Jxp', 0.)
        Jyp = model_params.get('Jyp', 0.)
        Jzp = model_params.get('Jzp', 0.)
        Jxpp = model_params.get('Jxpp', 0.)
        Jypp = model_params.get('Jypp', 0.)
        Jzpp = model_params.get('Jzpp', 0.)
        Jc = model_params.get('Jc', 0.)
        hx = model_params.get('hx', 0.)
        hy = model_params.get('hy', 0.)
        hz = model_params.get('hz', 0.)

        for u in range(len(self.lat.unit_cell)):
            self.add_onsite(-hx, u, 'Sx')
            self.add_onsite(-hy, u, 'Sy')
            self.add_onsite(-hz, u, 'Sz')
        for u1, u2, dx in self.lat.pairs['nearest_neighbors']:
            self.add_coupling((Jx + Jy) / 4., u1, 'Sp', u2, 'Sm', dx)
            # h.c.
            self.add_coupling(np.conj((Jx + Jy) / 4.), u2, 'Sp', u1, 'Sm', -dx)
            self.add_coupling((Jx - Jy) / 4., u1, 'Sp', u2, 'Sp', dx)
            # h.c.
            self.add_coupling(np.conj((Jx - Jy) / 4.), u2, 'Sm', u1, 'Sm', -dx)
            self.add_coupling(Jz, u1, 'Sz', u2, 'Sz', dx)
        for u1, u2, dx in self.lat.pairs['next_nearest_neighbors']:
            self.add_coupling((Jxp + Jyp) / 4., u1, 'Sp', u2, 'Sm', dx)
            # h.c.
            self.add_coupling(np.conj((Jxp + Jyp) / 4.),
                              u2, 'Sp', u1, 'Sm', -dx)
            self.add_coupling((Jxp - Jyp) / 4., u1, 'Sp', u2, 'Sp', dx)
            # h.c.
            self.add_coupling(np.conj((Jxp - Jyp) / 4.),
                              u2, 'Sm', u1, 'Sm', -dx)
            self.add_coupling(Jzp, u1, 'Sz', u2, 'Sz', dx)
        for u1, u2, dx in self.lat.pairs['next_next_nearest_neighbors']:
            self.add_coupling((Jxpp + Jypp) / 4., u1, 'Sp', u2, 'Sm', dx)
            # h.c.
            self.add_coupling(np.conj((Jxpp + Jypp) / 4.),
                              u2, 'Sp', u1, 'Sm', -dx)
            self.add_coupling((Jxpp - Jypp) / 4., u1, 'Sp', u2, 'Sp', dx)
            # h.c.
            self.add_coupling(np.conj((Jxpp - Jypp) / 4.),
                              u2, 'Sm', u1, 'Sm', -dx)
            self.add_coupling(Jzpp, u1, 'Sz', u2, 'Sz', dx)
        # Chiral interactions all clockwise
        for ii in range(len(Chiral_term)):
            # -i/2 S_z S_- S_+
            self.add_multi_coupling(-1j*Jc/2, [( 'Sz', [0, 0], Chiral_term[ii][0]), ( 'Sm', Chiral_term[ii][2], Chiral_term[ii][1]), ('Sp', Chiral_term[ii][4],Chiral_term[ii][3])])
            self.add_multi_coupling(np.conj(-1j*Jc/2), [( 'Sz', [0, 0], Chiral_term[ii][0]), ('Sp', Chiral_term[ii][2], Chiral_term[ii][1]), ( 'Sm', Chiral_term[ii][4], Chiral_term[ii][3])])
            # -i/2 S_+ S_z S_-
            self.add_multi_coupling(-1j*Jc/2, [( 'Sp', [0, 0], Chiral_term[ii][0]), ( 'Sz', Chiral_term[ii][2], Chiral_term[ii][1]), ('Sm', Chiral_term[ii][4],Chiral_term[ii][3])])
            self.add_multi_coupling(np.conj(-1j*Jc/2), [( 'Sm', [0, 0], Chiral_term[ii][0]), ('Sz', Chiral_term[ii][2], Chiral_term[ii][1]), ( 'Sp', Chiral_term[ii][4], Chiral_term[ii][3])])
            # -i/2 S_- S_+ S_z
            self.add_multi_coupling(-1j*Jc/2, [( 'Sm', [0, 0], Chiral_term[ii][0]), ( 'Sp', Chiral_term[ii][2], Chiral_term[ii][1]), ('Sz', Chiral_term[ii][4],Chiral_term[ii][3])])
            self.add_multi_coupling(np.conj(-1j*Jc/2), [( 'Sp', [0, 0], Chiral_term[ii][0]), ('Sm', Chiral_term[ii][2], Chiral_term[ii][1]), ( 'Sz', Chiral_term[ii][4], Chiral_term[ii][3])])

def DMRG_Square(Jxy, Jc, Lx, Ly, conserve='Sz', verbose=True):
    model_params = dict(
        lattice='Square',
        S=1,  # spin 1
        Lx=Lx,
        Ly=Ly,
        Jx=1,  # NN couplings
        Jy=1,
        Jz=1,
        Jxp=Jxy,  # nNN couplings
        Jyp=Jxy,
        Jzp=Jxy,
        Jxpp=0,  # nnNN couplings
        Jypp=0,
        Jzpp=0,
        Jc=Jc, # chiral interaction
        hx=0,
        hy=0,
        hz=0,
        bc_MPS='infinite',
        bc_y='cylinder',
        conserve=conserve,
        verbose=verbose)
    print(model_params['bc_MPS'], " DMRG, ", model_params['lattice'], " Spin-", model_params['S'])
    print("Jx={Jx:.4f}, Jxy={Jxy:.4f}, Jc={Jc:.4f}, conserve={conserve!r}".format(Jx=model_params['Jx'], 
        Jxy=Jxy, Jc=Jc, conserve=conserve))
    print("Lx={Lx:.2f}, Ly={Ly:.2f}".format(Lx=Lx, Ly=Ly))
    M = Chiral_term_for_Square(model_params)
    print('number of sites: ', M.lat.N_sites)

    # initial Neel state
    print("initial Neel state!", "\n")
    product_state = ["up", "down"] * int(model_params['Lx'] * model_params['Ly'] / 2)
    M.lat.test_sanity()
    psi = MPS.from_product_state(
        M.lat.mps_sites(), product_state, bc=M.lat.bc_MPS)

    dmrg_params = {
        'mixer': True,  # setting this to True helps to escape local minima
        'mixer_params': {
            'amplitude': 1.e-5,
            'decay': 1.2,
            'disable_after': 4
        },
        'trunc_params': {
            'svd_min': 1.e-10
            #'chi_max': 7000
        },
        'chi_list': {
            0: 80,
            4: 100,
            8: 200,
            12: 400,
            16: 800,
            20: 1500,
            24: 3000,
            28: 4000,
            32: 6000
        },
        'max_E_err': 1.e-9,
        'max_S_err': 1.e-6,
        'norm_tol': 1.e-6,
        'max_sweeps': 46,
        'max_hours': 324,
        'verbose': 1.,
        'N_sweeps_check': 2,
    }

    #DMRG starts
    print("DMRG starts", "\n")
    eng = dmrg.TwoSiteDMRGEngine(psi, M, dmrg_params)
    E, psi = eng.run()
    print("DMRG ends", "\n")
    print("E = {E:.13f}".format(E=E))
    print("final bond dimensions: ", psi.chi)
    Sz = psi.expectation_value("Sz")
    mag_z = np.mean(Sz)
    print("average ={mag_z:.9f}".format(mag_z=mag_z))
    stagger_mag_z = np.mean(np.absolute(Sz))
    print("absolute average ={stagger_mag_z:.9f}".format(
        stagger_mag_z=stagger_mag_z))
    if(model_params['bc_MPS'] == 'infinite'):
        xi = psi.correlation_length()
        print('correlation length =', xi, "\n")
    print("Site ; Spin ", "\n")
    for i in range(0, M.lat.N_sites):
        print(" ", (i+1), " ", Sz[i])
    print("\n")

    data = {}
    data['S'] = model_params['S']
    data['M'] = M
    data['E'] = E
    data['psi'] = psi
    data['bc_MPS'] = model_params['bc_MPS']
    data['bc_y'] = model_params['bc_y']
    data['lattice'] = model_params['lattice']
    data['Lx'] = Lx
    data['Ly'] = Ly
    data['Jxy'] = Jxy
    data['Jc'] = Jc

    return data

def Entanglement_spectrum_K_space(data):
    print("norm_test() = ", data['psi'].norm_test(), "\n")
    U, W, q, ov, trunc_err = data['psi'].compute_K(
        perm=data['M'].lat, verbose=5)
    print("Eigenvalue of TM = ", ov)
    print("truncation error = ", trunc_err)
    print("K space Entanglement spectrum: ")
    print("Sz ; eigenvalue ; K")
    for qi in range(q.block_number):
        K_spec = (q.get_charge(qi), np.sort(W[q.get_slice(qi)]))
        for i in range(0, len(K_spec[1])):
            if (np.abs(K_spec[1][i]) > 0.00003):
                print(K_spec[0][0], " ", np.abs(K_spec[1][i]),
                      " ", np.angle(K_spec[1][i]))
    print("\n")

if __name__ == "__main__":
    data = DMRG_Square(Jxy=0.45, Jc=0.3, Lx=1, Ly=10)  # N = Lx * (Ly)

    # K space entanglement spectrum
    if(data['bc_MPS'] == 'infinite'):
        Entanglement_spectrum_K_space(data)
    # Entanglement entropy at center bond
    EnEn = data['psi'].entanglement_entropy(n=1)
    print("Entanglement entropy at center= ",
          EnEn[int(data['M'].lat.N_sites / 2 - 1)], "\n")
    for i in range(0, len(EnEn)):
          print("Entanglement entropy at ", i, " = ", EnEn[i], "\n")
    
