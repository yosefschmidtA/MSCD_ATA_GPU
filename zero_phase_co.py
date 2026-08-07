import sys

def main():
    with open('psl9.txt', 'r') as f:
        lines = f.readlines()
        
    out_lines = []
    for line in lines:
        parts = line.split()
        if len(parts) == 10:
            try:
                # Check if the first is a float
                k = float(parts[0])
                # Format: k p(l=0) 0.0000 0.0000 ...
                new_parts = [parts[0], parts[1]] + ["0.0000"] * 8
                
                # Try to preserve spacing somewhat
                # Original format e.g.: "    5.0000   0.9957  -1.0914   2.6870..."
                out_lines.append("    {:<6}   {:<6}   {:<6}   {:<6}   {:<6}   {:<6}   {:<6}   {:<6}   {:<6}   {:<6}\n".format(*new_parts))
            except ValueError:
                out_lines.append(line)
        else:
            out_lines.append(line)
            
    with open('didatico_co.ph', 'w') as f:
        f.writelines(out_lines)

if __name__ == '__main__':
    main()
