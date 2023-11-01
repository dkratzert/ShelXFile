
class CifFile():
    """Class for writing IUCr CIF-1.1 of DDL1 compliant Version 2.4.5 files."""
    def __init__(self, shelx_file, template=None):
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

    def _write_cif(self):
        with open(self.template, "r") as f:
            template = f.read()
        cif = template.format(**self._cif_dict())
        return cif

    def _cif_dict(self):
        cif_dict = {}
        cif_dict["data_name"] = self.data.title
        cif_dict["cell"] = self._cell_dict()
        cif_dict["symmetry"] = self._symmetry_dict()
        cif_dict["atoms"] = self._atoms_dict()
        cif_dict["refine"] = self._refine_dict()
        cif_dict["special"] = self._special_dict()
        cif_dict["connections"] = self._connections_dict()
        cif_dict["misc"] = self._misc_dict()
        return cif_dict

    def _cell_dict(self):
        cell_dict = {}
        cell_dict["a"] = self.data.cell.a
        cell_dict["b"] = self.data.cell.b
        cell_dict["c"] = self.data.cell.c
        cell_dict["alpha"] = self.data.cell.alpha
        cell_dict["beta"] = self.data.cell.beta
        cell_dict["gamma"] = self.data.cell.gamma
        cell_dict["volume"] = self.data.cell.volume
        return cell_dict