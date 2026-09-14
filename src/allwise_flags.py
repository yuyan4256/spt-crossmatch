"""AllWISE per-band quality flags, decoded.

`qph` (photometric quality) and `ccf` (contamination/confusion) are 4-char
strings, one character per band W1/W2/W3/W4 (AllWISE Explanatory Supplement).
"""
BANDS = ('W1', 'W2', 'W3', 'W4')

QPH_MEANING = {
    'A': 'high SNR (≥10)', 'B': 'SNR 3-10', 'C': 'SNR 2-3',
    'U': 'upper limit (SNR<2)', 'X': 'no measurement',
    'Z': 'detected W1+W2 only',
}
# lower case = uncertain, upper case = certain
CCF_MEANING = {
    '0': 'no artifact', 'p': 'persistence (uncertain)', 'P': 'persistence',
    'd': 'diffraction spike (uncertain)', 'D': 'diffraction spike',
    'h': 'halo (uncertain)', 'H': 'halo',
    'o': 'optical ghost (uncertain)', 'O': 'optical ghost',
}


def decode(flags, table):
    """[(band, char, meaning)] for a 4-char flag string; [] if malformed."""
    if not isinstance(flags, str) or len(flags) != 4:
        return []
    return [(b, c, table.get(c, f'?{c}?')) for b, c in zip(BANDS, flags)]
