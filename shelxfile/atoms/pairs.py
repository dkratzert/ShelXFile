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
