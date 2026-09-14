"""Pure-numpy nested-HEALPix ang2pix (no healpy dependency).

Only what the Gaia tile downloader needs: ang2pix_nest for vectorized
(ra, dec) arrays. Algorithm follows the standard HEALPix primer /
healpy C implementation. Validated against Gaia DR3 source_id encoding
(source_id >> 35 == level-12 nested pixel of the source position).
"""
import numpy as np


def _spread_bits(v):
    """Insert a zero bit between each bit of v (v < 2^32)."""
    r = v.astype(np.uint64)
    r = (r | (r << np.uint64(16))) & np.uint64(0x0000FFFF0000FFFF)
    r = (r | (r << np.uint64(8))) & np.uint64(0x00FF00FF00FF00FF)
    r = (r | (r << np.uint64(4))) & np.uint64(0x0F0F0F0F0F0F0F0F)
    r = (r | (r << np.uint64(2))) & np.uint64(0x3333333333333333)
    r = (r | (r << np.uint64(1))) & np.uint64(0x5555555555555555)
    return r


def ang2pix_nest(nside, ra_deg, dec_deg):
    """Nested-scheme HEALPix pixel index for positions (degrees, ICRS)."""
    ra = np.atleast_1d(np.asarray(ra_deg, dtype=np.float64))
    dec = np.atleast_1d(np.asarray(dec_deg, dtype=np.float64))
    z = np.sin(np.deg2rad(dec))
    phi = np.deg2rad(ra % 360.0)
    za = np.abs(z)
    tt = (2.0 / np.pi) * phi % 4.0

    ix = np.empty(z.shape, dtype=np.int64)
    iy = np.empty(z.shape, dtype=np.int64)
    face = np.empty(z.shape, dtype=np.int64)

    eq = za <= 2.0 / 3.0
    if np.any(eq):
        temp1 = nside * (0.5 + tt[eq])
        temp2 = nside * (0.75 * z[eq])
        jp = np.floor(temp1 - temp2).astype(np.int64)
        jm = np.floor(temp1 + temp2).astype(np.int64)
        ifp = jp // nside
        ifm = jm // nside
        face[eq] = np.where(
            ifp == ifm, (ifp % 4) + 4, np.where(ifp < ifm, ifp % 4, (ifm % 4) + 8)
        )
        ix[eq] = jm & (nside - 1)
        iy[eq] = (nside - 1) - (jp & (nside - 1))

    po = ~eq
    if np.any(po):
        ntt = np.minimum(np.floor(tt[po]).astype(np.int64), 3)
        tp = tt[po] - ntt
        tmp = nside * np.sqrt(3.0 * (1.0 - za[po]))
        jp = np.floor(tp * tmp).astype(np.int64)
        jm = np.floor((1.0 - tp) * tmp).astype(np.int64)
        jp = np.minimum(jp, nside - 1)
        jm = np.minimum(jm, nside - 1)
        north = z[po] >= 0
        face[po] = np.where(north, ntt, ntt + 8)
        ix[po] = np.where(north, nside - jm - 1, jp)
        iy[po] = np.where(north, nside - jp - 1, jm)

    pix_in_face = (_spread_bits(ix) | (_spread_bits(iy) << np.uint64(1))).astype(
        np.int64
    )
    return face * nside * nside + pix_in_face
