def answer(l):
    for i in range(1, len(l) - 1):
        cnt_x = len([x for x in l[:i] if l[i] % x == 0])  # Count possible x in (x, y, z)
        cnt_z = len([z for z in l[i + 1:] if z % l[i] == 0])  # Count possible z in (x, y, z)

    return cnt_x * cnt_z


print answer([1, 1, 1])
print answer([1, 2, 3, 4, 5, 6])
