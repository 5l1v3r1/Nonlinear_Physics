import sys

if len(sys.argv) != 3:
    print ('Usage: $ chglast filename "line of text"')
    print()
    sys.exit()

filename = sys.argv[1]
line = sys.argv[2]

with open(filename, "a+") as f:
    f.seek(0, 2)
    pos = f.tell()
    print("Offset inital:", pos)
    pos -= 1
    f.seek(pos)
    # Dismiss end of line at end of file.
    s = f.read(1)
    if s == '\n':
        pos -= 1
    # Search previous end of line.
    f.seek(pos)
    while f.read(1) != '\n' and pos > 0:
        pos -= 1
        f.seek(pos)
    print("Offset final:", pos)
    f.seek(pos+1)
    f.truncate()
    f.write(line)
    if not line.endswith('\n'):
        f.write('\n')

print("done")

