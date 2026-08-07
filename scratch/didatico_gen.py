import os

input_file = '/home/yosef/MSCDATA/psAg111.txt'
output_file = '/home/yosef/MSCDATA/didatico.ph'

with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
    for i, line in enumerate(fin):
        # Keep header (first 12 lines typically, or until 'lnum')
        if i < 12:
            fout.write(line)
        else:
            parts = line.split()
            if not parts:
                fout.write(line)
                continue
            
            # Format the output line: k and l=0 remain, others become 0.0000
            try:
                k_val = float(parts[0])
                l0_val = float(parts[1])
                
                new_parts = [f"{k_val:10.4f}", f"{l0_val:9.4f}"]
                for _ in range(len(parts) - 2):
                    new_parts.append(f"{0.0:9.4f}")
                
                fout.write("  " + " ".join(new_parts) + "\n")
            except ValueError:
                # If not parsable, just write the line
                fout.write(line)

print("didatico.ph gerado com sucesso!")
