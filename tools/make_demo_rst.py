#!/usr/bin/env python3
"""Генератор синтетического несжатого ANSYS RST (little-endian, double).

Формат FORTRAN unformatted: [int32 nwords][body 4*nwords][int32 nwords].
Смещения в заголовках — в 4-байтных словах от начала файла.

Раскладка (совместима с парсером cpp/src/rst.cpp):
  rec@0      : PTRPTR, PTRDAT, PRSOL, PTRGEO, PTRELM, PTRCMS (f64->слова)
  @PTRDAT    : Glbnnod, _, Maxcomp
  @PRSOL     : NSETS, RPNTR(в словах), TPRNT3D(в словах)
  @RPNTR     : абсолютные смещения RSPEC каждого набора
  @TPRNT3D   : наборы по 3 int32 (step, substep, time)
  для набора s:
     RSPEC        : [_, DofBlkPtr(f64->слова)]
     DOF table    : тройки int32 (dtype, node, data_word); типы 9/10/11 = SX/SY/SZ
     data record  : 1 f64 (значение компонента для узла)
  @PTRCMS      : [NCMS,NNAM,NVAL][names][counts][node lists]
"""
import math, struct, sys

W = 4
def frame(body: bytes) -> bytes:
    n = len(body) // W
    assert len(body) % W == 0
    return struct.pack("<i", n) + body + struct.pack("<i", n)

class Builder:
    def __init__(self):
        self.buf = bytearray()
    @property
    def words(self): return len(self.buf) // W
    def add(self, body: bytes) -> int:
        off = self.words
        self.buf += frame(body)
        return off

def f64(*v): return struct.pack(f"<{len(v)}d", *v)
def i32(*v): return struct.pack(f"<{len(v)}i", *v)

N_NODES, N_SETS = 500, 8
NAME = b"FIXTURE".ljust(8)               # NNAM = 8
COMP_NODES = list(range(1, 101))

b = Builder()
rec1_off = b.add(f64(0, 0, 0, 0, 0, 0))          # placeholder, перезапишем
ptrdat = b.add(f64(N_NODES, 0, 1))

# --- компоненты ---
cms_body = i32(1, len(NAME), len(COMP_NODES)) + NAME + b"\x00" + i32(len(COMP_NODES)) + i32(*COMP_NODES)
if len(cms_body) % W: cms_body += b"\0" * (W - len(cms_body) % W)
ptrcms = b.add(cms_body)

# --- наборы результатов ---
rspec_offsets = []
for s in range(N_SETS):
    base = math.sin(0.1 * s)
    d_sx, d_sy, d_sz = [], [], []
    for nd in range(1, N_NODES + 1):
        d_sx.append(b.add(f64(base * nd)))
        d_sy.append(b.add(f64(base * nd * 0.5)))
        d_sz.append(b.add(f64(-base * nd * 0.25)))
    dof = bytearray()
    for nd in range(1, N_NODES + 1):
        dof += i32(9,  nd, d_sx[nd - 1])
        dof += i32(10, nd, d_sy[nd - 1])
        dof += i32(11, nd, d_sz[nd - 1])
    dof_off = b.add(bytes(dof))
    rspec_offsets.append(b.add(f64(0, dof_off)))

rpntr = b.add(f64(*rspec_offsets))
tprnt = bytearray()
for s in range(N_SETS):
    tprnt += i32(1, s + 1, s + 1)         # step, substep, time(int)
tprnt_off = b.add(bytes(tprnt))
prsolf = b.add(f64(N_SETS, rpntr, tprnt_off))

# --- перезапись первой записи: указатели таблиц ---
# Layout: [1]=PTRPTR, [2]=PTRDAT, [3]=PRSOL, [4]=PTRGEO, [5]=PTRELM, [6]=PTRCMS
new_rec1 = frame(f64(prsolf, ptrdat, prsolf, 0, 0, ptrcms))
start = rec1_off * W
old_len = len(frame(f64(0, 0, 0, 0, 0, 0)))
assert len(new_rec1) == old_len
b.buf[start:start + old_len] = new_rec1

out = sys.argv[1] if len(sys.argv) > 1 else "demo.rst"
with open(out, "wb") as fh:
    fh.write(bytes(b.buf))
print("wrote", out, len(b.buf), "bytes,", N_SETS, "sets,", N_NODES, "nodes")
