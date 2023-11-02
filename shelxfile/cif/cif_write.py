import datetime
from pathlib import Path
from string import Template
from typing import TYPE_CHECKING, Dict
from shelxfile import __version__
from shelxfile.atoms.atoms import Atoms

if TYPE_CHECKING:
    from shelxfile.shelx.shelxfile import Shelxfile


class CifFile():
    """Class for writing IUCr CIF-1.1 of DDL1 compliant Version 2.4.5 files."""

    def __init__(self, shelx_file: 'Shelxfile', template: str = None):
        if template is not None:
            self.template = template
        else:
            self.template = Path(__file__).parent / "cif_template.tmpl"
        self.data = shelx_file
        self._cif = None

    def __str__(self):
        return self._cif

    def __repr__(self):
        return self._cif

    def _write_cif(self) -> str:
        with open(self.template, "r") as f:
            template = f.read()
        cif = Template(template)
        sub = cif.substitute(self._cif_dict())
        return sub

    def _cif_dict(self) -> Dict[str, str]:
        cif_dict = {}
        cif_dict["data_name"] = self.data.titl.split()[0].lower()
        cif_dict["version"] = __version__
        cif_dict["creation_date"] = f"'{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}'"
        cif_dict["sum_formula"] = self.data.sum_formula
        cif_dict["formula_weight"] = 1234
        cif_dict.update(self._cell_data())
        cif_dict.update(self._symmetry_data())
        cif_dict.update(self._atoms_data())
        cif_dict.update(self._misc_dict())
        return cif_dict

    def _cell_data(self) -> Dict[str, float]:
        cell = self.data.cell
        return {
            "cell_a"     : cell.a,
            "cell_b"     : cell.b,
            "cell_c"     : cell.c,
            "cell_alpha" : cell.alpha,
            "cell_beta"  : cell.beta,
            "cell_gamma" : cell.gamma,
            "cell_volume": round(cell.volume, 4),
            "cell_z"     : self.data.zerr.Z,
        }

    def _symmetry_data(self) -> Dict[str, str]:
        symmcards_ = [f" '{sym.to_cif()}'" for sym in self.data.symmcards]
        return {
            "space_group"   : self.data.space_group,
            #"crystal_system": self.data.symmcards.crystal_system,
            "symmetry_loop" : '\n'.join(symmcards_),
        }

    def _atoms_data(self) -> Dict[str, str]:
        atoms: Atoms = self.data.atoms
        lines = []
        for atom in atoms.all_atoms:
            if not atom.qpeak:
                lines.append(f"{atom.name} {atom.element} "
                             f"{atom.x} {atom.y} {atom.z} "
                             f"{'Uiso' if atom.is_isotropic else 'Uani'} "
                             f"{atom.occupancy} {atom.part.n}")
        return {
            "atom_loop": '\n '.join(lines),
        }

    def _misc_dict(self) -> Dict[str, str]:
        misc_dict = {}
        misc_dict["temperature"] = round(self.data.temp_in_kelvin, 3)
        misc_dict["crystal_size_max"] = self.data.size.max
        misc_dict["crystal_size_mid"] = self.data.size.mid
        misc_dict["crystal_size_min"] = self.data.size.min
        misc_dict["wavelength"] = self.data.wavelength
        misc_dict["R1"] = self.data.R1
        misc_dict["wR2"] = self.data.wr2
        misc_dict["goodness_of_fit"] = self.data.goof
        return misc_dict


if __name__ == '__main__':
    from shelxfile import Shelxfile

    shx = Shelxfile(debug=True)
    shx.read_file('./tests/resources/p21c.res')
    # print(shx)
    cif = CifFile(shx)
    c = cif._write_cif()
    print(c)
    Path('p21c-test.cif').write_text(c)