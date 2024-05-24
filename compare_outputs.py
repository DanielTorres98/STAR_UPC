def compare_files(file1, file2):
    with open(file1, 'r') as f1, open(file2, 'r') as f2:
        lines1 = f1.readlines()
        lines2 = f2.readlines()

    max_lines = max(len(lines1), len(lines2))

    differences = []

    for i in range(max_lines):
        line1 = lines1[i].strip() if i < len(lines1) else "File1: <No Line>"
        line2 = lines2[i].strip() if i < len(lines2) else "File2: <No Line>"

        if line1 != line2:
            differences.append((i + 1, line1, line2))

    return differences

def main():
    file1 = 'output.txt'
    file2 = 'output_2.txt'

    differences = compare_files(file1, file2)

    if not differences:
        print("The files are identical.")
    else:
        print("Differences found:")
        for line_num, line1, line2 in differences:
            print(f"Line {line_num}:")
            print(f"File 1: {line1}")
            print(f"File 2: {line2}")
            print()

if __name__ == "__main__":
    main()
