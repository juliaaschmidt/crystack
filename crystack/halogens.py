

    # def is_functional_group(self, atom, group):
    #     """Given a pybel atom, look up if it belongs to a function group"""
    #     n_atoms = [a_neighbor.GetAtomicNum() for a_neighbor in 
    # pybel.ob.OBAtomAtomIter(atom.OBAtom)]

    #     if group in ['quartamine', 'tertamine'] and atom.atomicnum == 7:  # Nitrogen
    #         # It's a nitrogen, so could be a protonated amine or quaternary ammonium
    #         if '1' not in n_atoms and len(n_atoms) == 4:
    #             return True if group == 'quartamine' else False  # It's a quat. ammonium (N with 4 residues != H)
    #         elif atom.OBAtom.GetHyb() == 3 and len(n_atoms) >= 3:
    #             return True if group == 'tertamine' else False  # It's sp3-hybridized, so could pick up an hydrogen
    #         else:
    #             return False

    #     if group in ['sulfonium', 'sulfonicacid', 'sulfate'] and atom.atomicnum == 16:  # Sulfur
    #         if '1' not in n_atoms and len(n_atoms) == 3:  # It's a sulfonium (S with 3 residues != H)
    #             return True if group == 'sulfonium' else False
    #         elif n_atoms.count(8) == 3:  # It's a sulfonate or sulfonic acid
    #             return True if group == 'sulfonicacid' else False
    #         elif n_atoms.count(8) == 4:  # It's a sulfate
    #             return True if group == 'sulfate' else False

    #     if group == 'phosphate' and atom.atomicnum == 15:  # Phosphor
    #         if set(n_atoms) == {8}:  # It's a phosphate
    #             return True

    #     if group in ['carboxylate', 'guanidine'] and atom.atomicnum == 6:  # It's a carbon atom
    #         if n_atoms.count(8) == 2 and n_atoms.count(6) == 1:  # It's a carboxylate group
    #             return True if group == 'carboxylate' else False
    #         elif n_atoms.count(7) == 3 and len(n_atoms) == 3:  # It's a guanidine group
    #             nitro_partners = []
    #             for nitro in pybel.ob.OBAtomAtomIter(atom.OBAtom):
    #                 nitro_partners.append(len([b_neighbor for b_neighbor in pybel.ob.OBAtomAtomIter(nitro)]))
    #             if min(nitro_partners) == 1:  # One nitrogen is only connected to the carbon, can pick up a H
    #                 return True if group == 'guanidine' else False

    #     if group == 'halocarbon' and atom.atomicnum in [9, 17, 35, 53]:  # Halogen atoms
    #         n_atoms = [na for na in pybel.ob.OBAtomAtomIter(atom.OBAtom) if na.GetAtomicNum() == 6]
    #         if len(n_atoms) == 1:  # Halocarbon
    #             return True
    #     else:
    #         return False

    # def find_hal_acc(self, atoms):
    #     """Look for halogen bond acceptors (Y-{O|P|N|S}, with Y=C,P,S)"""
    #     data = namedtuple('hal_acceptor', 'o o_orig_idx y y_orig_idx')
    #     a_set = []
    #     # All oxygens, nitrogen, sulfurs with neighboring [carbon, phosphor, nitrogen or sulfur]
    #     for a in [at for at in atoms if at.atomicnum in [8, 7, 16]]:
    #         n_atoms = [na for na in pybel.ob.OBAtomAtomIter(a.OBAtom) if na.GetAtomicNum() in [6, 7, 15, 16]]
    #         if len(n_atoms) == 1:  # Proximal atom
    #             #o_orig_idx = self.Mapper.mapid(a.idx, mtype=self.mtype, bsid=self.bsid)
    #             #y_orig_idx = self.Mapper.mapid(n_atoms[0].GetIdx(), mtype=self.mtype, bsid=self.bsid)
    #             a_set.append(data(o=a, o_orig_idx=None, y=pybel.Atom(n_atoms[0]), y_orig_idx=None))
    #     return a_set

    # def find_hal_donors(self, atoms):
    #     """Look for halogen bond donors (X-C, with X=F, Cl, Br, I)"""
    #     data = namedtuple('hal_donor', 'x orig_x x_orig_idx c c_orig_idx')
    #     a_set = []
    #     for a in atoms:
    #         if Interactions.is_functional_group(self, a, 'halocarbon'):
    #             n_atoms = [na for na in pybel.ob.OBAtomAtomIter(a.OBAtom) if na.GetAtomicNum() == 6]
    #             #x_orig_idx = self.Mapper.mapid(a.idx, mtype=self.mtype, bsid=self.bsid)
    #             #orig_x = self.Mapper.id_to_atom(x_orig_idx)
    #             #c_orig_idx = [self.Mapper.mapid(na.GetIdx(), mtype=self.mtype, bsid=self.bsid) for na in n_atoms]
    #             a_set.append(data(x=a, orig_x=None, x_orig_idx=None,
    #                               c=pybel.Atom(n_atoms[0]), c_orig_idx=None))
    #     # if len(a_set) != 0:
    #     #     write_message('Ligand contains %i halogen atom(s).\n' % len(a_set), indent=True)
    #     return a_set

    # def halogen(acceptor, donor):
    #     """Detect all halogen bonds of the type Y-O...X-C"""
    #     data = namedtuple('halogenbond', 'acc acc_orig_idx don don_orig_idx distance don_angle acc_angle restype '
    #                                      'resnr reschain restype_l resnr_l reschain_l donortype acctype sidechain')
    #     pairings = []
    #     for acc, don in itertools.product(acceptor, donor):
    #         dist = euclidean3d(acc.o.coords, don.x.coords)
    #         if not config.MIN_DIST < dist < config.HALOGEN_DIST_MAX:
    #             continue
    #         vec1, vec2 = vector(acc.o.coords, acc.y.coords), vector(acc.o.coords, don.x.coords)
    #         vec3, vec4 = vector(don.x.coords, acc.o.coords), vector(don.x.coords, don.c.coords)
    #         acc_angle, don_angle = vecangle(vec1, vec2), vecangle(vec3, vec4)
    #         is_sidechain_hal = acc.o.OBAtom.GetResidue().GetAtomProperty(acc.o.OBAtom, 8)  # Check if sidechain atom
    #         if not config.HALOGEN_ACC_ANGLE - config.HALOGEN_ANGLE_DEV < acc_angle \
    #                 < config.HALOGEN_ACC_ANGLE + config.HALOGEN_ANGLE_DEV:
    #             continue
    #         if not config.HALOGEN_DON_ANGLE - config.HALOGEN_ANGLE_DEV < don_angle \
    #                 < config.HALOGEN_DON_ANGLE + config.HALOGEN_ANGLE_DEV:
    #             continue

    #         contact = data(acc=acc, acc_orig_idx=acc.o_orig_idx, don=don, don_orig_idx=don.x_orig_idx,
    #                        distance=dist, don_angle=don_angle, acc_angle=acc_angle,
    #                        restype=None, resnr=None,
    #                        reschain=None, restype_l=None,
    #                        reschain_l=None, resnr_l=None, donortype=don.x.OBAtom.GetType(), acctype=acc.o.type,
    #                        sidechain=None)
    #         print(f"Halogen interaction found: {don.x.OBAtom.GetType()} {acc.o.type} {dist}")
    #         pairings.append(contact)
    #     return pairings #filter_contacts(pairings)


    def calculate_halogen_interactions(self):
        #currently directional feature: oxygens on one molecule, halogens on other molecule
        halogenbond_acc1 = Interactions.find_hal_acc(self, self.MolA)
        halogenbond_don1 = Interactions.find_hal_donors(self, self.MolA)

        halogenbond_acc2 = Interactions.find_hal_acc(self, self.MolB)
        halogenbond_don2 = Interactions.find_hal_donors(self, self.MolB)

        # Thus consider both directions: oxygens on one molecule, halogens on other molecule AND vice versa
        halogen_bonds1 = Interactions.halogen(halogenbond_acc1, halogenbond_don2)
        halogen_bonds2 = Interactions.halogen(halogenbond_acc2, halogenbond_don1)

        halogen_bonds = [halogen_bonds1, halogen_bonds2]