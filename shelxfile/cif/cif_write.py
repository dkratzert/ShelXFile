import datetime
from pathlib import Path
from string import Template
from typing import TYPE_CHECKING
from shelxfile import __version__

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

    def _get_cif(self):
        if self._cif is None:
            self._cif = self._write_cif()
        return self._cif

    def _write_cif(self) -> Template:
        with open(self.template, "r") as f:
            template = f.read()
        cif = Template(template)
        sub = cif.substitute(self._cif_dict())
        return sub

    def _cif_dict(self):
        cif_dict = {}
        cif_dict["data_name"] = self.data.titl.split()[0].lower()
        cif_dict["version"] = __version__
        cif_dict["creation_date"] = f"'{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}'"
        cif_dict["sum_formula"] = self.data.sum_formula
        cif_dict["formula_weight"] = 1234
        cif_dict.update(self._cell_data())
        cif_dict.update(self._symmetry_data())
        cif_dict["atoms"] = self._atoms_data()
        # cif_dict["refine"] = self._refine_data()
        # cif_dict["misc"] = self._misc_dict()
        return cif_dict

    def _cell_data(self):
        cell = self.data.cell
        return {
            "cell_a"     : cell.a,
            "cell_b"     : cell.b,
            "cell_c"     : cell.c,
            "cell_alpha" : cell.alpha,
            "cell_beta"  : cell.beta,
            "cell_gamma" : cell.gamma,
            "cell_volume": round(cell.volume, 4),
        }

    def _symmetry_data(self):
        print(self.data.symmcards.latt_ops)
        return {
            "space_group"  : self.data.space_group,
            "crystal_system": self.data.symmcards.latt_ops,
            # "hall_symbol"              : symmetry.hall_symbol,
            # "origin"                   : symmetry.origin,
            # "centering"                : symmetry.centering,
            # "unit_cell_setting"        : symmetry.unit_cell_setting,
            "symmetry_loop": self.data.symmcards.latt_ops,
        }

    def _atoms_data(self):
        atoms = self.data.atoms
        return {
            "atoms": atoms,
        }


if __name__ == '__main__':
    from shelxfile import Shelxfile

    shx = Shelxfile(debug=True)
    shx.read_file('./tests/resources/p21c.res')
    # print(shx)
    cif = CifFile(shx)
    c = cif._write_cif()
    print(c)
