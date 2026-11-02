"""Minimal Renishaw WDF reader: returns (x, spectra) from XLST and DATA blocks."""
import struct
import numpy as np

def read_wdf(path):
    b = open(path, "rb").read()
    assert b[:4] == b"WDF1", "not a WDF file"
    # header: 512 bytes; then blocks (4s tag, uint32 uid, uint64 size incl. 16-byte header)
    npoints = struct.unpack_from("<I", b, 60)[0]      # points per spectrum
    ncoll = struct.unpack_from("<Q", b, 64)[0]        # number of spectra
    off, blocks = 512, {}
    while off + 16 <= len(b):
        tag, uid, size = struct.unpack_from("<4sIQ", b, off)
        if size < 16: break
        blocks[tag.decode(errors="replace")] = (off + 16, size - 16)
        off += size
    xo, xs = blocks["XLST"]
    x = np.frombuffer(b, "<f4", count=npoints, offset=xo + 8)  # skip uint32 type + uint32 unit
    do, ds = blocks["DATA"]
    y = np.frombuffer(b, "<f4", count=npoints * ncoll, offset=do).reshape(ncoll, npoints)
    return x.astype(float), y.astype(float), blocks
