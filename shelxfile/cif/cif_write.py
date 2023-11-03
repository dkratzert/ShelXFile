import datetime
from pathlib import Path
from string import Template
from typing import TYPE_CHECKING, Dict
from shelxfile.version import VERSION
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

    def write_cif(self, cif_path: Path) -> None:
        with open(self.template, "r") as f:
            template = f.read()
        cif = Template(template)
        sub = cif.substitute(self._cif_dict())
        cif_path.write_text(sub)

    def _cif_dict(self) -> Dict[str, str]:
        cif_dict = {}
        cif_dict["data_name"] = self.data.titl.split()[0].lower() or "unknown"
        cif_dict["version"] = VERSION
        cif_dict["creation_date"] = f"{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}"
        cif_dict["sum_formula"] = self.data.sum_formula
        cif_dict["formula_weight"] = self.data.formula_weight
        cif_dict.update(self._cell_data())
        cif_dict.update(self._symmetry_data())
        cif_dict.update(self._atoms_data())
        cif_dict.update(self._adp_data())
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
            "space_group"  : self.data.space_group,
            # "crystal_system": self.data.symmcards.crystal_system,
            "symmetry_loop": '\n'.join(symmcards_),
        }

    def _atoms_data(self) -> Dict[str, str]:
        atoms: Atoms = self.data.atoms
        lines = []
        loop_header = ("loop_\n"
                       " _atom_site_label\n"
                       " _atom_site_type_symbol\n"
                       " _atom_site_fract_x\n"
                       " _atom_site_fract_y\n"
                       " _atom_site_fract_z\n"
                       " _atom_site_adp_type\n"
                       " _atom_site_occupancy\n"
                       " _atom_site_disorder_group")
        for atom in atoms.all_atoms:
            if not atom.qpeak:
                lines.append(f"{atom.fullname_short} {atom.element} "
                             f"{atom.x} {atom.y} {atom.z} "
                             f"{'Uiso' if atom.is_isotropic else 'Uani'} "
                             f"{atom.occupancy} {atom.part.n}")
        return {
            "atom_loop_header": loop_header,
            "atom_loop"       : '\n'.join(lines),
        }

    def _adp_data(self) -> Dict[str, str]:
        atoms: Atoms = self.data.atoms
        lines = []
        loop_header = (
            "loop_\n"
            " _atom_site_aniso_label\n"
            " _atom_site_aniso_U_11\n"
            " _atom_site_aniso_U_22\n"
            " _atom_site_aniso_U_33\n"
            " _atom_site_aniso_U_23\n"
            " _atom_site_aniso_U_13\n"
            " _atom_site_aniso_U_12")
        for atom in atoms.all_atoms:
            if not atom.qpeak and not atom.is_isotropic:
                lines.append(f"{atom.fullname_short} "
                             f"{atom.U11} {atom.U22} {atom.U33} "
                             f"{atom.U23} {atom.U13} {atom.U12} ")
        return {
            "aniso_loop_header": loop_header,
            "aniso_loop"       : '\n'.join(lines),
        }

    def _misc_dict(self) -> Dict[str, str]:
        misc_dict = {}
        misc_dict["temperature"] = round(self.data.temp_in_kelvin, 3) or '?'
        misc_dict["crystal_size_max"] = self.data.size.max or '?'
        misc_dict["crystal_size_mid"] = self.data.size.mid or '?'
        misc_dict["crystal_size_min"] = self.data.size.min or '?'
        misc_dict["wavelength"] = self.data.wavelength or '?'
        misc_dict["R1"] = self.data.R1 or '?'
        misc_dict["wR2"] = self.data.wr2 or '?'
        misc_dict["goodness_of_fit"] = self.data.goof or '?'
        return misc_dict


if __name__ == '__main__':
    from shelxfile import Shelxfile

    shx = Shelxfile(debug=True)
    #shx.read_file('tests/resources/I-43d.res')
    shx.read_file('./tests/resources/p21c.res')

    #cif = CifFile(shx)
    #cif.write_cif(Path('p21c-test.cif'))
    shx.to_cif('p21c-test.cif')
