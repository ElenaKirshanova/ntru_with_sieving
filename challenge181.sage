q = 1024
N = 181

h = [955, 826, 662, 668, 264, 927, 771, 509, 929, 97, 533,
401, 845, 896, 53, 426, 127, 413, 253, 614, 948, 376, 273, 264, 901,
701, 539, 331, 225, 390, 206, 668, 287, 977, 100, 185, 683, 638, 792,
904, 944, 311, 984, 558, 404, 983, 269, 362, 812, 88, 484, 760, 827,
820, 554, 881, 900, 290, 315, 821, 735, 881, 374, 542, 679, 1020, 562,
661, 233, 239, 872, 465, 646, 507, 221, 125, 416, 78, 32, 180, 819,
946, 753, 961, 1015, 326, 3, 774, 564, 404, 785, 786, 394, 624, 258,
90, 256, 578, 302, 580, 610, 245, 200, 165, 213, 225, 975, 799, 302,
393, 116, 717, 942, 566, 384, 240, 623, 556, 600, 869, 771, 225, 692,
1019, 949, 8, 693, 400, 379, 698, 77, 369, 22, 644, 955, 828, 174,
640, 211, 119, 720, 926, 907, 32, 832, 232, 1018, 106, 65, 841, 273,
295, 921, 984, 64, 147, 373, 990, 761, 922, 628, 581, 994, 594, 718,
529, 600, 125, 352, 953, 231, 274, 557, 151, 827, 309, 835, 654, 884,
339, 514]
sol = [-1, 0, -1, -1, -1, -1, 0, 1, -1, 1, 1, 0, -1, -1, 0, 1, 0, 1, 1, -1, 1, 1, 0, 1, -1, 1, 1, 0, 0, 0, -1, 1, 0, -1, 0, -1, -1, -1, 0, -1, 0, -1, 0, 1, 1, 1, 1, -1, -1, 1, -1, 0, -1, -1, 0, 1, 1, 1, 1, -1, -1, 0, 0, 0, -1, -1, 0, 0, -1, 1,
1, 0, 1, 0, -1, 0, 1, 0, 0, 1, 1, -1, -1, 1, -1, 1, 0, 0, 0, -1, 0, 0, -1, -1, 1, 1, -1, 0, 0, 0, 1, -1, 1, 0, -1, 1, 1, -1, -1, 0, -1, -1, 0, -1, 0, 1, -1, 1, 0, 0, 1, -1, 0, 1, 0, -1, -1, 0, 1, 0, 1, -1, 1, -1, -1, -1, 0, 1, 1, -1, 1, 0, 0, -1, -1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, -1, 1, 1, 1, -1, -1, 0, -1, 1, 0, 0, 0, 1, 1, 1, -1, -1, 1, -1, 0, -1, 1, 1, 1, 0, 1, 1, 0, -1, -1, 0, 0, -1, 0, 0, -1, -1, 0, -1, -1, -1, 1, -1, 0, -1, -1, -2, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, -1, 2, 0, 0, 1, 2, 1, 0, -1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 2, 2, 1, 0, 0, 0, 1, 2, 3, 0, 0, 1, 0, 0, 0, 0, 1, -2, -2, 0, -1, 0, 0, 0, 0, -2, 1, 1, -1, 1, 0, 0, 0, -1, -1, 0, -1, -2, -1, 0, 0, 0,
0, 0, -1, -1, 0, 1, 0, 0, 0, 1, -1, -1, -1, 2, -2, -1, 0, 0, -1, -1, -1, 1, -2, 1, -1, 0, 0, 0, 1, 1, -1, 0, 0, 0, 1, 0, 0, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 1, 2, 2, 0, 1, 1, -1, 0, 0, 2, 1, 1, 2, -1, 1, 0, 0, 1, 0, -1, 1, 0, 0, 0, 0, 0, 0, -1]

g = sol[:181]
f = sol[181:]
Zq     = Integers(q)
Zqx    = PolynomialRing(Zq, 'x')
Zqxquo = Zqx.quotient(x^N - 1, 'x')

gpoly = Zqxquo(g)
fpoly = Zqxquo(f)
hpoly = Zqxquo(h)

#print(sum(g))
print(gpoly)
print(-(1-3*fpoly)*hpoly)
