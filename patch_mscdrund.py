import sys

with open("mscdrund_cpu.cpp", "r") as f:
    orig = f.read()

with open("mscdrund.cpp", "r") as f:
    curr = f.read()

orig_start = orig.find("  MSCDT_AMARK;\n  for (ib=0;(error==0)&&(ib<natoms);++ib)")
orig_end = orig.find("  MSCDT_A(4,\"summation: laco de m (serie)\");") + len("  MSCDT_A(4,\"summation: laco de m (serie)\");\n")
cpu_code = orig[orig_start:orig_end]

curr_start = curr.find("#ifdef MSCDGPU\n  MSCDT_AMARK;\n  error = mscdgpu_summation")
curr_end = curr.find("  MSCDT_A(4,\"summation: laco de m (gpu)\");\n#else\n") + len("  MSCDT_A(4,\"summation: laco de m (gpu)\");\n#else\n")

if curr_start != -1:
    # We replace the whole block we added before
    curr_end_full = curr.find("  MSCDT_A(4,\"summation: laco de m (serie)\");\n#endif\n") + len("  MSCDT_A(4,\"summation: laco de m (serie)\");\n#endif\n")
    
    new_block = """#ifdef MSCDGPU
  if (getenv("MSCD_GPU")) {
    MSCDT_AMARK;
    error = mscdgpu_summation(akin, (const Gcplx*)tevenelem, (Gcplx*)asum, patom);
    MSCDT_A(4,"summation: laco de m (gpu)");
  } else {
#endif
""" + cpu_code + """#ifdef MSCDGPU
  }
#endif
"""
    new_curr = curr[:curr_start] + new_block + curr[curr_end_full:]
    with open("mscdrund.cpp", "w") as f:
        f.write(new_curr)
    print("Patched!")
else:
    print("Could not find block in mscdrund.cpp")
