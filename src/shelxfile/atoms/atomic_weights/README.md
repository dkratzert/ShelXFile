# atomicWeightsDecimal

A Python dictionary of atomic weights in Decimal objects.

Measurements come from IUPAC's Comission on Isotopic Abundances and Atomic Weights (IUPAC-CIAAW). This dictionary uses measurements from IUPAC-CIAWW's **Atomic weights of the elements: Review 2000**.

See [ciaaw.org](https://www.caaw.og) for updated measurements.

## Implementation

Per IUPAC-CIAAW, only elements with stable isotopes are included.

Using an element's chemical symbol (e.g., `"H"`) as a key returns a sub-dictionary. The sub-dictionary contains two keys (1) `standard` for the atomic weight and (2) `abundant` for the atomic mass of the element's most-abundant isotope.

`Standard atomic weight` or `atomic weight` refers to the weighted average mass of an element's isotopes. 

The `most-abundant isoptope` of an element is important for high resolution mass sepctromentry. For biologically relevant molecules, an elements most-abundant isotope is usually the lightest isotope. An `exact mass` can be calculated from these values.

Tab characters are used as whitespace in the dictionary to simplify regex manipulation.

A for loop iterates through all chemical element's weights.

## Caution

Care should be taken when applying these measurements. Standard weights may not be applicable to high-resolution mass spectrometry. Some elements have multiple stable isotopes with similar relative abundances. For instance, 50.69% of Bromine is <sup>79</sup>Br and the other 49.31% is <sup>81</sup>Br. Chlorine's isotope abundances are particularily relevant to biochemists: 75.76% <sup>35</sup>Cl (`34.9688527 71 Da`) and 24.24% <sup>37</sup>Cl (`36.96590260 Da`).