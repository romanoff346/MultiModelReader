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
# Элементная named selection: quad4-элементы из узлов 201..300 (сетка 10x10)
ELEM_NAME = b"ELEMS".ljust(8)

b = Builder()
rec1_off = b.add(f64(0, 0, 0, 0, 0, 0, 0))       # placeholder (7 слотов), перезапишем
ptrdat = b.add(f64(N_NODES, 0, 2))               # Maxcomp = 2

# --- геометрия/узлы: @PTRGEO [NELMAX][NELEM]... (упрощённо — только счётчики) ---
# Таблица элементов: ELMAP = [NELMAX][NELEM], ELCONN = id + 8 нод на элемент
elems = []            # (eid, [n1..n4])
nid = lambda r, c: 201 + r * 10 + c
for r in range(9):
    for c in range(9):
        elems.append((len(elems) + 1, [nid(r, c), nid(r, c + 1),
                                       nid(r + 1, c + 1), nid(r + 1, c)]))
elconn_body = bytearray()
for eid, conn in elems:
    row = conn + [0] * (8 - len(conn))
    elconn_body += i32(*row)                      # 8 слов нод на элемент
elmap_off_tmp = b.words                           # адрес ELMAP известен заранее
elmap_body = i32(len(elems), len(elems)) + f64(elmap_off_tmp + 1)  # ptr ELCONN (в словах)
elmap_off = b.add(elmap_body)
elconn_off = b.add(bytes(elconn_body))
ptrelm = b.add(f64(elmap_off))                    # PTRELM -> таблица указателей эл.

# --- компоненты (узловой + элементный) ---
def comp_block(names, types, lists):
    ncms, nnam = len(names), max(len(n) for n in names)
    body = i32(ncms, nnam, sum(len(l) for l in lists))
    for nm in names: body += nm.ljust(nnam) + b"\x00"
    for t in types:  body += t.ljust(4)           # маркеры типов 'N'/'E'
    pad = (-len(body)) % W                        # выравнивание на слово
    body += b"\x00" * pad
    for l in lists:  body += i32(len(l)) + i32(*l)
    if len(body) % W: body += b"\x00" * (W - len(body) % W)
    return body

cms_body = comp_block([NAME, ELEM_NAME], [b"N", b"E"], [COMP_NODES, [e[0] for e in elems]])
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
new_rec1 = frame(f64(prsolf, ptrdat, prsolf, 0, ptrelm, ptrcms, elmap_off))
start = rec1_off * W
old_len = len(frame(f64(0, 0, 0, 0, 0, 0, 0)))
assert len(new_rec1) == old_len
b.buf[start:start + old_len] = new_rec1

out = sys.argv[1] if len(sys.argv) > 1 else "demo.rst"
with open(out, "wb") as fh:
    fh.write(bytes(b.buf))
print("wrote", out, len(b.buf), "bytes,", N_SETS, "sets,", N_NODES, "nodes")
