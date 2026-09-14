class AtomPair():
    def __init__(self, atom1, atom2):
        self.atom1 = atom1
        self.atom2 = atom2

    def __repr__(self):
        if self.atom1 and self.atom2:
            return f'{self.atom1} {self.atom2}'
        else:
            return ''

    def __len__(self):
        if self.atom1 and self.atom2:
            return 2
        elif not self.atom1 or not self.atom2:
            return 0


class Bond:
    """A covalent bond between two atoms with a measured distance.

    Attributes
    ----------
    atom1, atom2 : :class:`~shelxfile.atoms.atom.Atom`
        The bonded atom pair (atom1 always has the lower ``all_atoms`` index).
    distance : float
        Interatomic distance in Å.
    """

    def __init__(self, atom1, atom2, distance: float) -> None:
        self.atom1 = atom1
        self.atom2 = atom2
        self.distance = distance

    def __iter__(self):
        """Allows unpacking: ``a1, a2, dist = bond``."""
        yield self.atom1
        yield self.atom2
        yield self.distance

    def __repr__(self) -> str:
        return (f"Bond({self.atom1.fullname_short}, "
                f"{self.atom2.fullname_short}, "
                f"{self.distance:.4f} Å)")

    def __str__(self) -> str:
        return (f"{self.atom1.fullname_short:<6s}"
                f" – "
                f"{self.atom2.fullname_short:<6s}"
                f"  {self.distance:.4f} Å")


class SymBond(Bond):
    """A covalent bond that may involve a symmetry-generated neighbor.

    Extends :class:`Bond` with a *symm_label* that describes the symmetry
    operation (and any lattice translation) needed to generate the neighbor
    atom from its position in the asymmetric unit.  An empty label means
    both atoms are in the asymmetric unit with no additional translation.

    The label uses SHELXL / CIF convention, e.g. ``'-x, y+1/2, -z+1/2'``.

    Attributes
    ----------
    atom1 : :class:`~shelxfile.atoms.atom.Atom`
        The asymmetric-unit atom that "owns" the bond (the central atom).
    atom2 : :class:`~shelxfile.atoms.atom.Atom`
        The neighbor in the asymmetric unit (its image is generated via
        *symm_label* to obtain the actual bonded position).
    distance : float
        Interatomic distance in Å.
    symm_label : str
        Symmetry operation string, e.g. ``'-x, y+1/2, -z+1/2'``, or an
        empty string for plain asymmetric-unit bonds.
    symm_number : int
        0-based index of the symmetry operation in ``Shelxfile.symmcards``
        (0 = identity).  Always 0 for plain asymmetric-unit bonds.
    """

    def __init__(self, atom1, atom2, distance: float,
                 symm_label: str = '', symm_number: int = 0) -> None:
        super().__init__(atom1, atom2, distance)
        self.symm_label = symm_label
        self.symm_number = symm_number

    @property
    def is_symmetry_bond(self) -> bool:
        """``True`` when the bond involves a symmetry-generated image."""
        return bool(self.symm_label)

    def __iter__(self):
        """Allows unpacking: ``a1, a2, dist, label = bond``."""
        yield self.atom1
        yield self.atom2
        yield self.distance
        yield self.symm_label

    def __repr__(self) -> str:
        if self.symm_label:
            return (f"SymBond({self.atom1.fullname_short}, "
                    f"{self.atom2.fullname_short}, "
                    f"{self.distance:.4f} Å, "
                    f"#{self.symm_number} [{self.symm_label}])")
        return (f"SymBond({self.atom1.fullname_short}, "
                f"{self.atom2.fullname_short}, "
                f"{self.distance:.4f} Å)")

    def __str__(self) -> str:
        suffix = ''
        if self.symm_label:
            suffix = f'   #{self.symm_number} [{self.symm_label}]'
        return (f"{self.atom1.fullname_short:<6s}"
                f" – "
                f"{self.atom2.fullname_short:<6s}"
                f"  {self.distance:.4f} Å"
                f"{suffix}")
